//! Functions to apply on meta cells. These are in-memory versions of the
//! variants used for single cells due to the massive compression achieved
//! with meta cells.

use bixverse_rs::prelude::*;
use bixverse_rs::single_cell::mc_analysis::metrics::pairwise_gene_correlations_in_memory;
use bixverse_rs::single_cell::mc_processing::hvg_pca::*;
use bixverse_rs::single_cell::sc_data::in_memory_io::InMemorySparseReader;
use bixverse_rs::single_cell::sc_processing::analytic_pearson::model::AprParams;
use bixverse_rs::single_cell::sc_processing::analytic_pearson::stream::{
    apr_gene_pass, build_apr_model, cell_totals_over_genes, fit_analytic_pearson_grouped,
    AprGroupedFit,
};
use bixverse_rs::single_cell::sc_processing::hvg::*;
use bixverse_rs::single_cell::sc_processing::pca::{pca_on_sc_residuals, SingleCellPcaParams};
use bixverse_rs::single_cell::sc_processing::residuals::residual_variance;
use bixverse_rs::single_cell::sc_processing::sctransform::model::{SctCovariates, SctParams};
use bixverse_rs::single_cell::sc_processing::sctransform::stream::{
    fit_sctransform, fit_sctransform_grouped, SctGroupedFit, SctStreamOpts,
};
use extendr_api::*;
use std::time::Instant;

use crate::meta_cell::utils::{mc_list_to_sparse_f32, mc_list_to_sparse_u32};
use crate::single_cell::r_sc_residuals::{
    apr_fit_to_r_list, parse_residual_fit, residual_variance_to_r_list, resolve_groups,
    sct_fit_to_r_list, stream_opts, with_residual_source, METHOD_APR, METHOD_SCT,
};

/////////////
// extendR //
/////////////

extendr_module! {
    mod r_mc_processing;
    // hvg and pca
    fn rs_mc_hvg;
    fn rs_mc_pca;
    // residuals
    fn rs_mc_fit_residuals;
    fn rs_mc_residual_variance;
    fn rs_mc_pca_residuals;
    // correlations
    fn rs_pairwise_gene_cors_mc;
}

///////////////////////////
// Highly variable genes //
///////////////////////////

/// Meta cells highly variable genes
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Calculates highly variable genes for MetaCells or more
/// generally speaking sparse data. This is happening in-memory compared to the
/// (usually much) larger single cell data sets. `"meanvarbin"` and
/// `"dispersion"` compute the same statistics; they differ only in how the R
/// side selects from them.
///
/// @param sparse_data A named list that needs to have `data`, `indptr`,
/// `indices`, `nrow`, `ncol` and `cs_type`. Shape is (metacells, genes). Pass
/// raw counts for `"vst"` and normalised counts otherwise.
/// @param hvg_method String. Which HVG detection method to use. Options
/// are `c("vst", "meanvarbin", "dispersion")`.
/// @param loess_span Numeric. The span parameter for the loess function
/// (only used for `"vst"`).
/// @param binning String. The binning strategy for the `meanvarbin` and
/// `dispersion` methods. One of `c("equal_width", "equal_freq")`.
/// @param n_bins Integer. Number of bins for the `meanvarbin` and
/// `dispersion` methods.
/// @param clip_max Optional clipping number. Defaults to `sqrt(no_cells)` if
/// not provided (only used for `"vst"`).
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @returns A list with the HVG statistics. If `hvg_method == "vst"`:
/// \itemize{
///   \item mean - The average expression of the gene.
///   \item var - The variance of the gene.
///   \item var_exp - The expected variance of the gene.
///   \item var_std - The standardised variance of the gene.
/// }
/// For `"meanvarbin"` and `"dispersion"`:
/// \itemize{
///   \item mean - The average expression of the gene.
///   \item dispersion - The dispersion of the gene.
///   \item dispersion_scaled - The scaled dispersion per bin per gene.
///   \item bin - The bin of the gene.
/// }
///
/// @export
#[extendr]
fn rs_mc_hvg(
    sparse_data: List,
    hvg_method: &str,
    loess_span: f64,
    binning: String,
    n_bins: usize,
    clip_max: Option<f32>,
    verbose: usize,
) -> Result<List> {
    let start = Instant::now();
    let verbosity = parse_verbosity_level(verbose);

    if verbosity.normal_verbosity() {
        println!("Running HVG detection on meta cells.")
    }

    let sparse = mc_list_to_sparse_f32(sparse_data)?;

    if verbosity.detailed_verbosity() {
        println!(
            " MetaCell HVG: Finished data transformation in {:.2?}...",
            start.elapsed()
        )
    }

    let hvg_type = parse_hvg_method(hvg_method)
        .ok_or_else(|| format!("Invalid HVG method: {}", hvg_method))
        .unwrap();

    let res = match hvg_type {
        HvgMethod::Vst => {
            let res = get_hvg_vst_from_sparse(&sparse, loess_span as f32, clip_max);
            Ok(list!(
                mean = res.mean,
                var = res.var,
                var_exp = res.var_exp,
                var_std = res.var_std
            ))
        }
        HvgMethod::MeanVarBin => {
            let res = get_hvg_mvb_from_sparse(&sparse, &binning, n_bins);
            Ok(list!(
                mean = res.mean,
                dispersion = res.dispersion,
                dispersion_scaled = res.dispersion_scaled,
                bin = res.bin
            ))
        }
        HvgMethod::Dispersion => {
            let res = get_hvg_dispersion_from_sparse(&sparse, &binning, n_bins);
            Ok(list!(
                mean = res.mean,
                dispersion = res.dispersion,
                dispersion_scaled = res.dispersion_scaled,
                bin = res.bin
            ))
        }
    };

    if verbosity.detailed_verbosity() {
        println!(
            " MetaCell HVG: Finished HVG calculations in {:.2?}...",
            start.elapsed()
        )
    }

    res
}

/////////
// PCA //
/////////

/// PCA on MetaCells (sparse data)
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Calculates PCA for MetaCells or more generally speaking sparse
/// data. This is happening in-memory compared to the (usually much) larger
/// single cell data sets. The matrix is densified, optionally CLR transformed
/// and scaled according to `pca_params` before the SVD.
///
/// @param sparse_data A named list that needs to have `data`, `indptr`,
/// `indices`, `nrow`, `ncol` and `cs_type`. Shape is (metacells, genes),
/// holding the normalised counts of the genes to use.
/// @param no_pcs Integer. Number of PCs to return.
/// @param pca_params Named list. Contains the parameters to use for this PCA
/// run, see [bixverse::params_sc_pca()].
/// @param clr_offsets Optional numeric. One offset per meta cell for the
/// `PFlogPF` normalisation from Booeshaghi, et al., computed against the full
/// gene panel. Required if `pca_params$clr` is `TRUE`, ignored otherwise.
/// @param seed Integer. Random seed for the randomised SVD.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @returns A list with the following items
/// \itemize{
///   \item scores - The samples projected on the PCA space (solved via sparse
///   SVD).
///   \item loadings - The loadings of the features for the PCA (solved via
///   sparse SVD).
///   \item singular_values - The singular values for the PCA (solved via sparse
///   SVD).
/// }
///
/// @export
///
/// @references Booeshaghi, et al., bioRxive, 2026.
#[extendr]
fn rs_mc_pca(
    sparse_data: List,
    no_pcs: usize,
    pca_params: List,
    clr_offsets: Option<Vec<f64>>,
    seed: usize,
    verbose: usize,
) -> Result<List> {
    let start = Instant::now();
    let verbosity = parse_verbosity_level(verbose);

    if verbosity.normal_verbosity() {
        println!("Running PCA calculation on meta cells.")
    }

    let sparse = mc_list_to_sparse_f32(sparse_data)?;

    if verbosity.detailed_verbosity() {
        println!(
            " MetaCell PCA: Finished data transformation in {:.2?}...",
            start.elapsed()
        )
    }

    let pca_params = SingleCellPcaParams::from_r_list(pca_params)?;

    let offsets = if pca_params.clr {
        let offsets = clr_offsets.ok_or_else(|| Error::Other("'clr_offsets' ".into()))?;
        Some(offsets)
    } else {
        None
    };

    let res =
        pca_on_metacells(&sparse, no_pcs, &pca_params, offsets.as_deref(), seed).to_extendr()?;

    if verbosity.detailed_verbosity() {
        println!(
            " MetaCell PCA: Finished SVD calculation in {:.2?}...",
            start.elapsed()
        )
    }

    Ok(list!(
        scores = faer_to_r_matrix(res.0.as_ref()),
        loadings = faer_to_r_matrix(res.1.as_ref()),
        singular_values = res.2.r_float_convert()
    ))
}

///////////////
// Residuals //
///////////////

/// Fits a residual model for meta cells
///
/// @description
/// `r lifecycle::badge("experimental")`
/// In-memory version of [bixverse::rs_sc_fit_residuals()]. Meta cell counts are
/// summed UMIs, so the negative binomial the residual models describe is still
/// defined; it just sits at a much greater depth than a single cell.
///
/// @param sparse_data A named list that needs to have `data`, `indptr`,
/// `indices`, `nrow`, `ncol` and `cs_type`. Shape is (metacells, genes). Pass
/// raw counts.
/// @param method String. One of `c("sctransform", "analytic_pearson")`.
/// @param group_of_cell Integer vector or `NULL`. Group label per meta cell.
/// (0-indexed, dense!) `NULL` fits one model.
/// @param covariates Named list of numeric vectors, one per covariate, each of
/// length `nrow`. scTransform only.
/// @param params Named list. See [bixverse::params_sc_sctransform()] or
/// [bixverse::params_sc_apr()].
/// @param seed Integer. Seed for the step-1 subsample.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @returns A list as described in [bixverse::rs_sc_fit_residuals()].
///
/// @export
///
/// @keywords internal
#[allow(clippy::too_many_arguments)]
#[extendr]
fn rs_mc_fit_residuals(
    sparse_data: List,
    method: &str,
    group_of_cell: Option<Vec<i32>>,
    covariates: List,
    params: List,
    seed: u64,
    verbose: usize,
) -> Result<List> {
    let sparse = mc_list_to_sparse_u32(sparse_data)?;
    let sparse = to_csc(sparse);

    let n_cells = sparse.shape.0;
    let cell_set: Vec<usize> = (0..n_cells).collect();
    let opts = stream_opts(None, seed, verbose);

    let gene_reader = InMemorySparseReader::new(&sparse, None).to_extendr()?;

    match method {
        METHOD_SCT => {
            let sct_params = SctParams::from_r_list(params)?;
            let covariates = SctCovariates::from_r_list(covariates)?;
            covariates.validate(n_cells).to_extendr()?;

            let library_sizes: Vec<f64> = gene_reader
                .read_cell_library_sizes(&cell_set)
                .to_extendr()?
                .into_iter()
                .map(|s| s as f64)
                .collect();
            let log10_umi: Vec<f64> = library_sizes.iter().map(|&s| s.log10()).collect();

            let fit = match &group_of_cell {
                None => {
                    let (model, _pass) = fit_sctransform(
                        &gene_reader,
                        &cell_set,
                        &library_sizes,
                        &covariates,
                        &sct_params,
                        None,
                        None,
                        opts,
                    )
                    .to_extendr()?;

                    SctGroupedFit {
                        genes: model.genes.clone(),
                        models: vec![model],
                        passes: Vec::new(),
                        group_of_cell: vec![0_u32; n_cells],
                    }
                }
                Some(_) => {
                    let groups = resolve_groups(group_of_cell.clone(), n_cells)?;
                    fit_sctransform_grouped(
                        &gene_reader,
                        &cell_set,
                        &groups,
                        &library_sizes,
                        &covariates,
                        &sct_params,
                        opts,
                    )
                    .to_extendr()?
                }
            };

            Ok(sct_fit_to_r_list(&fit, &log10_umi, &covariates, &cell_set))
        }
        METHOD_APR => {
            let apr_params = AprParams::from_r_list(params)?;

            // the analytic Pearson fit sums each cell over its group's retained
            // genes, which a gene-major sweep cannot answer
            let cell_reader =
                InMemorySparseReader::new_cell_major(&sparse, None).to_extendr()?;

            let fit = match &group_of_cell {
                None => {
                    let pass =
                        apr_gene_pass(&gene_reader, &cell_set, &apr_params, opts).to_extendr()?;
                    let totals = cell_totals_over_genes(&cell_reader, &cell_set, &pass.retained)
                        .to_extendr()?;
                    let model = build_apr_model(&pass, &totals, &apr_params).to_extendr()?;

                    AprGroupedFit {
                        models: vec![model],
                        cell_totals: totals,
                        group_of_cell: vec![0_u32; n_cells],
                    }
                }
                Some(_) => {
                    let groups = resolve_groups(group_of_cell.clone(), n_cells)?;
                    fit_analytic_pearson_grouped(
                        &gene_reader,
                        &cell_reader,
                        &cell_set,
                        &groups,
                        &apr_params,
                        opts,
                    )
                    .to_extendr()?
                }
            };

            Ok(apr_fit_to_r_list(&fit, &cell_set))
        }
        other => Err(Error::Other(format!(
            "Unknown residual method '{other}'. Expected '{METHOD_SCT}' or '{METHOD_APR}'."
        ))),
    }
}

/// Residual variance and variable features for meta cells
///
/// @description
/// `r lifecycle::badge("experimental")`
/// In-memory version of [bixverse::rs_sc_residual_variance()].
///
/// @param sparse_data A named list that needs to have `data`, `indptr`,
/// `indices`, `nrow`, `ncol` and `cs_type`. Shape is (metacells, genes). Pass
/// raw counts.
/// @param residual_fit List. A fit from [bixverse::rs_mc_fit_residuals()].
/// @param n_hvg Integer. Variable features to take from each group.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @returns A list with the following items
/// \itemize{
///   \item genes - The gene indices the variances are indexed by. (0-indexed!)
///   \item variance - Matrix of residual variance, genes by groups.
///   \item hvg - The selected gene indices, ascending. (0-indexed!)
/// }
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_mc_residual_variance(
    sparse_data: List,
    residual_fit: List,
    n_hvg: usize,
    verbose: usize,
) -> Result<List> {
    let sparse = mc_list_to_sparse_u32(sparse_data)?;
    let sparse = to_csc(sparse);

    let cell_set: Vec<usize> = (0..sparse.shape.0).collect();
    let opts = stream_opts(None, SctStreamOpts::default().seed, verbose);

    let reader = InMemorySparseReader::new(&sparse, None).to_extendr()?;
    let fit = parse_residual_fit(residual_fit, &cell_set)?;

    let (genes, variance, hvg) = with_residual_source(&fit, |source| {
        let per_group = residual_variance(&reader, source, &cell_set, opts)?;
        let hvg = select_residual_hvg(&per_group, source.genes(), n_hvg)?;
        Ok((source.genes().to_vec(), per_group, hvg))
    })?;

    Ok(residual_variance_to_r_list(&genes, &variance, &hvg))
}

/// Calculates PCA on Pearson residuals for meta cells
///
/// @description
/// `r lifecycle::badge("experimental")`
/// In-memory version of [bixverse::rs_sc_pca_residuals()]. As there, `clr` and
/// `normalise_variance` must both be `FALSE`.
///
/// @param sparse_data A named list that needs to have `data`, `indptr`,
/// `indices`, `nrow`, `ncol` and `cs_type`. Shape is (metacells, genes). Pass
/// raw counts.
/// @param residual_fit List. A fit from [bixverse::rs_mc_fit_residuals()].
/// @param no_pcs Integer. Number of PCs to calculate.
/// @param pca_params Named list. Contains the parameters to use for this PCA
/// run.
/// @param gene_indices Integer vector. The gene indices to use. (0-indexed!)
/// @param seed Integer. Random seed for the randomised SVD.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @returns A list with the following items
/// \itemize{
///   \item scores - The samples projected on the PCA space.
///   \item loadings - The loadings of the features for the PCA.
///   \item singular_values - The singular values for the PCA.
/// }
///
/// @export
///
/// @keywords internal
#[allow(clippy::too_many_arguments)]
#[extendr]
fn rs_mc_pca_residuals(
    sparse_data: List,
    residual_fit: List,
    no_pcs: usize,
    pca_params: List,
    gene_indices: Vec<i32>,
    seed: usize,
    verbose: usize,
) -> Result<List> {
    let sparse = mc_list_to_sparse_u32(sparse_data)?;
    let sparse = to_csc(sparse);

    let cell_set: Vec<usize> = (0..sparse.shape.0).collect();
    let gene_set = gene_indices.r_int_convert();

    let pca_params = SingleCellPcaParams::from_r_list(pca_params)?;

    let reader = InMemorySparseReader::new(&sparse, None).to_extendr()?;
    let fit = parse_residual_fit(residual_fit, &cell_set)?;

    let res = with_residual_source(&fit, |source| {
        pca_on_sc_residuals(
            &reader,
            &cell_set,
            &gene_set,
            no_pcs,
            &pca_params,
            source,
            seed,
            false,
            verbose,
        )
    })?;

    Ok(list!(
        scores = faer_to_r_matrix(res.0.as_ref()),
        loadings = faer_to_r_matrix(res.1.as_ref()),
        singular_values = res.2.r_float_convert()
    ))
}

/// Ensure the meta cell matrix is CSC over (metacells, genes)
///
/// The R side hands over CSR, and [InMemorySparseReader] wants the gene-major
/// layout of the same matrix. `transform` converts the format and leaves the
/// shape alone, which is exactly that; `transpose_and_convert` would swap the
/// axes instead.
///
/// ### Params
///
/// * `sparse` - The matrix as it arrived from R.
///
/// ### Returns
///
/// The same matrix, CSC.
fn to_csc(sparse: CompressedSparseData2<u32, f32>) -> CompressedSparseData2<u32, f32> {
    if sparse.cs_type.is_csc() {
        sparse
    } else {
        sparse.transform()
    }
}

///////////////////
// Pairwise cors //
///////////////////

/// Calculate the pairwise gene-correlation for meta cells
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Correlates `gene_indices_1[i]` against `gene_indices_2[i]` over the meta
/// cells, in memory.
///
/// @param sparse_data A named list that needs to have `data`, `indptr`,
/// `indices`, `nrow`, `ncol` and `cs_type`. Shape is (metacells, genes),
/// holding the normalised counts.
/// @param gene_indices_1 Integer. The gene indices for the first set of genes.
/// Must be 0-indexed!
/// @param gene_indices_2 Integer. The gene indices for the second set of
/// genes, same length as `gene_indices_1`. Must be 0-indexed!
/// @param spearman Boolean. Shall the Spearman correlation be calculated.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @returns Numeric vector with one correlation per pair of `gene_indices_1`
/// and `gene_indices_2`.
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_pairwise_gene_cors_mc(
    sparse_data: List,
    gene_indices_1: &[i32],
    gene_indices_2: &[i32],
    spearman: bool,
    verbose: usize,
) -> Result<Vec<f64>> {
    let gene_indices_1 = gene_indices_1.r_int_convert();
    let gene_indices_2 = gene_indices_2.r_int_convert();

    let sparse = mc_list_to_sparse_f32(sparse_data)?;

    // transpose
    let sparse = if sparse.cs_type.is_csr() {
        sparse.transform()
    } else {
        sparse
    };

    let pairwise_cors = pairwise_gene_correlations_in_memory(
        &sparse,
        &gene_indices_1,
        &gene_indices_2,
        spearman,
        verbose,
    )
    .to_extendr()?;

    Ok(pairwise_cors.r_float_convert())
}
