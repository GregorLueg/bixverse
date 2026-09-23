//! Rust <> R interface for the residual-based single cell path.
//!
//! scTransform and analytic Pearson both fit a null model and then regenerate
//! residual rows on demand, so the fitted model is the thing that travels
//! between calls. It goes to R as a plain named list and comes back the same
//! way, which keeps the R side stateless about how the residuals are computed
//! and lets one fit serve the HVG selection, the PCA and the count correction.
//!
//! The two methods differ only in the fit. Everything downstream is written
//! against `ResidualSource`, so the variance, PCA and correction bindings take
//! whichever model the list carries.

use bixverse_rs::prelude::*;
use bixverse_rs::single_cell::sc_data::data_io::ParallelSparseReader;
use bixverse_rs::single_cell::sc_processing::analytic_pearson::model::{AprParams, AprResiduals};
use bixverse_rs::single_cell::sc_processing::analytic_pearson::stream::{
    apr_gene_pass, build_apr_model, cell_totals_over_genes, fit_analytic_pearson_grouped,
    AprGroupedFit,
};
use bixverse_rs::single_cell::sc_processing::hvg::select_residual_hvg;
use bixverse_rs::single_cell::sc_processing::pca::{pca_on_sc_residuals, SingleCellPcaParams};
use bixverse_rs::single_cell::sc_processing::residuals::{
    intersect_gene_sets, residual_variance, ResidualSource,
};
use bixverse_rs::single_cell::sc_processing::sctransform::model::{
    SctCellContext, SctCovariates, SctParams,
};
use bixverse_rs::single_cell::sc_processing::sctransform::residuals::SctResiduals;
use bixverse_rs::single_cell::sc_processing::sctransform::stream::{
    fit_sctransform, fit_sctransform_grouped, sct_corrected_counts, SctGroupedFit, SctStreamOpts,
};
use bixverse_rs::single_cell::sc_data::bin_merge_io::gene_store_to_cell_store;
use extendr_api::prelude::*;

/////////////
// extendR //
/////////////

extendr_module! {
    mod r_sc_residuals;
    // fitting
    fn rs_sc_fit_residuals;
    // hvg
    fn rs_sc_residual_variance;
    // pca
    fn rs_sc_pca_residuals;
    // corrected counts and the transpose that makes them usable
    fn rs_sct_corrected_counts;
    fn rs_sc_gene_store_to_cell_store;
}

///////////////
// Constants //
///////////////

/// Name the R list uses for the scTransform method.
pub(crate) const METHOD_SCT: &str = "sctransform";

/// Name the R list uses for the analytic Pearson method.
pub(crate) const METHOD_APR: &str = "analytic_pearson";

///////////
// Types //
///////////

/// A fitted residual model, parsed back out of the R list.
///
/// Owns everything the matching [`ResidualSource`] borrows. The source itself
/// cannot be stored next to its inputs, so this is the owned half and
/// [`with_residual_source`] builds the borrowing half on the stack.
pub(crate) enum RResidualFit {
    /// scTransform, one model per group.
    Sct {
        /// The per-group models plus the shared gene axis.
        fit: SctGroupedFit,
        /// `log10(library size)` per selected cell, the model's offset.
        log10_umi: Vec<f64>,
        /// Cell-level covariates, in the order the models were fitted with.
        covariates: SctCovariates,
    },
    /// Analytic Pearson, one model per group.
    Apr {
        /// The per-group models plus the per-cell totals.
        fit: AprGroupedFit,
    },
}

//////////////
// Fit <> R //
//////////////

/// Serialise a fitted scTransform model set for R.
///
/// ### Params
///
/// * `fit` - The grouped fit, one model per group.
/// * `log10_umi` - Offset per selected cell.
/// * `covariates` - The covariates the models were fitted with.
/// * `cell_indices` - The selected cells, echoed so a later call can check it
///   is being handed the same selection.
///
/// ### Returns
///
/// The list [`parse_residual_fit`] reads back.
pub(crate) fn sct_fit_to_r_list(
    fit: &SctGroupedFit,
    log10_umi: &[f64],
    covariates: &SctCovariates,
    cell_indices: &[usize],
) -> List {
    let mut out = fit.to_r_list();
    let echoed: Vec<i32> = cell_indices.iter().map(|&c| c as i32).collect();

    out = extend_list(
        out,
        &[
            ("method", Robj::from(METHOD_SCT)),
            ("cell_indices", Robj::from(echoed)),
            ("covariates", Robj::from(covariates.to_r_list())),
            ("log10_umi", Robj::from(log10_umi.to_vec())),
        ],
    );

    out
}

/// Serialise a fitted analytic Pearson model set for R.
///
/// ### Params
///
/// * `fit` - The grouped fit, one model per group.
/// * `cell_indices` - The selected cells, echoed for the same reason as above.
///
/// ### Returns
///
/// The list [`parse_residual_fit`] reads back.
pub(crate) fn apr_fit_to_r_list(fit: &AprGroupedFit, cell_indices: &[usize]) -> List {
    let echoed: Vec<i32> = cell_indices.iter().map(|&c| c as i32).collect();

    // The analytic Pearson fit derives its shared gene axis in
    // `AprResiduals::new`, but R wants it for the HVG bookkeeping, so compute
    // it once here rather than making every caller intersect the models.
    let sets: Vec<&[usize]> = fit.models.iter().map(|m| m.genes.as_slice()).collect();
    let genes: Vec<i32> = intersect_gene_sets(&sets)
        .map(|g| g.iter().map(|&x| x as i32).collect())
        .unwrap_or_default();

    extend_list(
        fit.to_r_list(),
        &[
            ("method", Robj::from(METHOD_APR)),
            ("cell_indices", Robj::from(echoed)),
            ("genes", Robj::from(genes)),
        ],
    )
}

/// Append named elements to a list.
///
/// `list!()` cannot take a runtime-built set of names, and the grouped fits
/// already produce most of the list, so the R-side extras are appended rather
/// than the whole thing being rebuilt by hand.
///
/// ### Params
///
/// * `base` - The list to extend.
/// * `extra` - Name and value pairs to append.
///
/// ### Returns
///
/// The extended list.
fn extend_list(base: List, extra: &[(&str, Robj)]) -> List {
    let mut names: Vec<String> = base
        .names()
        .map(|n| n.map(String::from).collect())
        .unwrap_or_default();
    let mut values: Vec<Robj> = base.values().collect();

    for (name, value) in extra {
        names.push((*name).to_string());
        values.push(value.clone());
    }

    let mut out = List::from_values(values);
    out.set_names(names.iter().map(|s| s.as_str()))
        .expect("one name per value by construction");
    out
}

/// Rebuild a fitted model from the list the fit binding returned.
///
/// Also checks the fit belongs to the selection it is about to be used on. A
/// cell filter that moved after fitting would otherwise give plausible
/// residuals computed against the wrong cells, which nothing downstream can
/// detect.
///
/// ### Params
///
/// * `fit` - The serialised fit.
/// * `cell_indices` - The cells the caller is about to use it on.
///
/// ### Returns
///
/// The owned fit, or an error naming what does not line up.
pub(crate) fn parse_residual_fit(fit: List, cell_indices: &[usize]) -> Result<RResidualFit, extendr_api::Error> {
    let method: String = fit
        .index("method")
        .map_err(|_| Error::Other("The residual fit is missing 'method'".to_string()))?
        .as_str()
        .ok_or_else(|| Error::Other("The residual fit 'method' is not a string".to_string()))?
        .to_string();

    let fitted_cells: Vec<usize> = fit
        .index("cell_indices")
        .map_err(|_| Error::Other("The residual fit is missing 'cell_indices'".to_string()))?
        .as_integer_slice()
        .ok_or_else(|| {
            Error::Other("The residual fit 'cell_indices' is not integer".to_string())
        })?
        .iter()
        .map(|&c| c as usize)
        .collect();

    if fitted_cells != cell_indices {
        return Err(Error::Other(format!(
            "The residual fit was fitted on {} cells that are not the {} cells \
             given here. Refit with fit_residuals_sc() after changing the cell \
             selection.",
            fitted_cells.len(),
            cell_indices.len()
        )));
    }

    match method.as_str() {
        METHOD_SCT => {
            let log10_umi = fit
                .index("log10_umi")
                .map_err(|_| {
                    Error::Other("The scTransform fit is missing 'log10_umi'".to_string())
                })?
                .as_real_vector()
                .ok_or_else(|| {
                    Error::Other("The scTransform fit 'log10_umi' is not numeric".to_string())
                })?;

            let covariates = match fit.index("covariates") {
                Ok(robj) if !robj.is_null() => {
                    let list = robj.as_list().ok_or_else(|| {
                        Error::Other("The scTransform fit 'covariates' is not a list".to_string())
                    })?;
                    SctCovariates::from_r_list(list)?
                }
                _ => SctCovariates::default(),
            };

            Ok(RResidualFit::Sct {
                fit: SctGroupedFit::from_r_list(fit)?,
                log10_umi,
                covariates,
            })
        }
        METHOD_APR => Ok(RResidualFit::Apr {
            fit: AprGroupedFit::from_r_list(fit)?,
        }),
        other => Err(Error::Other(format!(
            "Unknown residual method '{other}'. Expected '{METHOD_SCT}' or '{METHOD_APR}'."
        ))),
    }
}

/// Run a closure against the residual source a fit describes.
///
/// Both `SctCellContext` and `SctResiduals` borrow their inputs, so a function
/// handing back a `&dyn ResidualSource` cannot be written: the owned data has
/// to outlive the borrow, which means it lives in this frame and the caller's
/// work happens inside the closure.
///
/// ### Params
///
/// * `fit` - The owned fit.
/// * `f` - What to do with the source.
///
/// ### Returns
///
/// Whatever the closure returned, or the first error either side raised.
pub(crate) fn with_residual_source<R>(
    fit: &RResidualFit,
    f: impl FnOnce(&dyn ResidualSource) -> std::result::Result<R, BixverseErrors>,
) -> Result<R, extendr_api::Error> {
    match fit {
        RResidualFit::Sct {
            fit,
            log10_umi,
            covariates,
        } => {
            let cells = SctCellContext::new(log10_umi, covariates).to_extendr()?;
            let source =
                SctResiduals::new(&fit.models, cells, fit.group_of_cell.clone()).to_extendr()?;
            f(&source).to_extendr()
        }
        RResidualFit::Apr { fit } => {
            let source =
                AprResiduals::new(&fit.models, &fit.cell_totals, fit.group_of_cell.clone())
                    .to_extendr()?;
            f(&source).to_extendr()
        }
    }
}

/// Stream options from the arguments every binding takes.
///
/// ### Params
///
/// * `gene_batch_size` - Genes held in memory per batch, `None` for the crate
///   default.
/// * `seed` - Seed for the step-1 subsample.
/// * `verbose` - `0` quiet, `1` normal, `2` detailed.
///
/// ### Returns
///
/// The [SctStreamOpts].
pub(crate) fn stream_opts(gene_batch_size: Option<usize>, seed: u64, verbose: usize) -> SctStreamOpts {
    let defaults = SctStreamOpts::default();

    SctStreamOpts {
        gene_batch_size: gene_batch_size.or(defaults.gene_batch_size),
        seed,
        verbose,
    }
}

/// Group labels for a fit, defaulting every cell into one group.
///
/// ### Params
///
/// * `group_of_cell` - The 0-based labels from R, or `None`.
/// * `n_cells` - Number of selected cells.
///
/// ### Returns
///
/// One label per selected cell, or an error when the length disagrees.
pub(crate) fn resolve_groups(group_of_cell: Option<Vec<i32>>, n_cells: usize) -> Result<Vec<u32>, extendr_api::Error> {
    match group_of_cell {
        None => Ok(vec![0_u32; n_cells]),
        Some(groups) => {
            if groups.len() != n_cells {
                return Err(Error::Other(format!(
                    "'group_of_cell' has {} entries but {n_cells} cells were selected",
                    groups.len()
                )));
            }
            groups
                .into_iter()
                .map(|g| {
                    u32::try_from(g).map_err(|_| {
                        Error::Other("'group_of_cell' has a negative entry".to_string())
                    })
                })
                .collect()
        }
    }
}

/////////////
// Fitting //
/////////////

/// Fits a residual model for single cell data
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Fits either scTransform (v2) or the analytic Pearson residual model of
/// Lause, Berens and Kobak over the selected cells. With `group_of_cell` one
/// model is fitted per group, which is what a multi-sample experiment wants:
/// each sample keeps its own sequencing depth and composition.
///
/// The fit is returned as a list rather than applied to anything. Hand it back
/// to [bixverse::rs_sc_residual_variance()], [bixverse::rs_sc_pca_residuals()]
/// or [bixverse::rs_sct_corrected_counts()] to use it.
///
/// @param f_path_gene String. Path to the `counts_genes.bin` file.
/// @param f_path_cell String. Path to the `counts_cells.bin` file. Used for the
/// library sizes, and for the per-cell totals of the analytic Pearson fit.
/// @param method String. One of `c("sctransform", "analytic_pearson")`.
/// @param cell_indices Integer vector. The cell indices to use. (0-indexed!)
/// @param group_of_cell Integer vector or `NULL`. Group label per selected
/// cell. (0-indexed, dense!) `NULL` fits one model over every cell.
/// @param covariates Named list of numeric vectors, one per covariate, each of
/// length `length(cell_indices)`. scTransform only. The order is remembered and
/// checked on every subsequent use.
/// @param params Named list. The parameters, see
/// [bixverse::params_sc_sctransform()] or [bixverse::params_sc_apr()].
/// @param gene_batch_size Integer or `NULL`. Genes held in memory per batch.
/// @param seed Integer. Seed for the step-1 subsample.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @returns A list with the following items
/// \itemize{
///   \item method - The method that was fitted.
///   \item models - The per-group models.
///   \item genes - The gene indices modelled in every group. (0-indexed!)
///   \item group_of_cell - The group label per selected cell. (0-indexed!)
///   \item n_groups - Number of groups.
///   \item cell_indices - The cells that were fitted on. (0-indexed!)
///   \item covariates - The covariates used, scTransform only.
///   \item log10_umi - The per-cell offset, scTransform only.
///   \item cell_totals - The per-cell totals, analytic Pearson only.
/// }
///
/// @export
///
/// @references Choudhary and Satija, Genome Biology, 2022; Lause, Berens and
/// Kobak, Genome Biology, 2021.
///
/// @keywords internal
#[allow(clippy::too_many_arguments)]
#[extendr]
fn rs_sc_fit_residuals(
    f_path_gene: &str,
    f_path_cell: &str,
    method: &str,
    cell_indices: Vec<i32>,
    group_of_cell: Option<Vec<i32>>,
    covariates: List,
    params: List,
    gene_batch_size: Option<usize>,
    seed: u64,
    verbose: usize,
) -> Result<List, extendr_api::Error> {
    let cell_set = cell_indices.r_int_convert();
    let n_cells = cell_set.len();
    let opts = stream_opts(gene_batch_size, seed, verbose);

    let gene_reader = ParallelSparseReader::new(f_path_gene).to_extendr()?;
    let cell_reader = ParallelSparseReader::new(f_path_cell).to_extendr()?;

    match method {
        METHOD_SCT => {
            let sct_params = SctParams::from_r_list(params)?;
            let covariates = SctCovariates::from_r_list(covariates)?;
            covariates.validate(n_cells).to_extendr()?;

            // The library sizes are a property of the store, so they are read
            // here rather than asked of R, where a stale vector could silently
            // shift every offset.
            let library_sizes: Vec<f64> = cell_reader
                .read_cell_library_sizes(&cell_set)
                .to_extendr()?
                .into_iter()
                .map(|s| s as f64)
                .collect();
            let log10_umi: Vec<f64> = library_sizes.iter().map(|&s| s.log10()).collect();

            let fit = match &group_of_cell {
                // A single fit is not the one-group case of the grouped fit:
                // it subsamples for step 1 over the whole selection.
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

            let fit = match &group_of_cell {
                None => {
                    let pass =
                        apr_gene_pass(&gene_reader, &cell_set, &apr_params, opts).to_extendr()?;
                    let totals =
                        cell_totals_over_genes(&cell_reader, &cell_set, &pass.retained)
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

/////////
// HVG //
/////////

/// Residual variance and the variable features it selects
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Computes the per-gene residual variance from a fitted model and selects the
/// variable features from it. With more than one group the selection follows
/// Seurat: rank within each group, take the top `n_hvg` of each, and union
/// them, so a marker only one sample carries is not buried by a pooled ranking.
/// The returned set can therefore be larger than `n_hvg`.
///
/// Each gene's residual row is regenerated, reduced and dropped, so memory is
/// one row per worker rather than a genes-by-cells matrix.
///
/// @param f_path_gene String. Path to the `counts_genes.bin` file.
/// @param residual_fit List. A fit from [bixverse::rs_sc_fit_residuals()].
/// @param cell_indices Integer vector. The cell indices to use. (0-indexed!)
/// Must be the selection the fit was fitted on.
/// @param n_hvg Integer. Variable features to take from each group.
/// @param gene_batch_size Integer or `NULL`. Genes held in memory per batch.
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
/// @references Seurat v5, `SCTransform.StdAssay`
///
/// @keywords internal
#[extendr]
fn rs_sc_residual_variance(
    f_path_gene: &str,
    residual_fit: List,
    cell_indices: Vec<i32>,
    n_hvg: usize,
    gene_batch_size: Option<usize>,
    verbose: usize,
) -> Result<List, extendr_api::Error> {
    let cell_set = cell_indices.r_int_convert();
    let opts = stream_opts(gene_batch_size, SctStreamOpts::default().seed, verbose);

    let reader = ParallelSparseReader::new(f_path_gene).to_extendr()?;
    let fit = parse_residual_fit(residual_fit, &cell_set)?;

    let (genes, variance, hvg) = with_residual_source(&fit, |source| {
        let per_group = residual_variance(&reader, source, &cell_set, opts)?;
        let hvg = select_residual_hvg(&per_group, source.genes(), n_hvg)?;
        Ok((source.genes().to_vec(), per_group, hvg))
    })?;

    Ok(residual_variance_to_r_list(&genes, &variance, &hvg))
}

/// Pack the residual variance sweep into the list R expects.
///
/// ### Params
///
/// * `genes` - The gene axis the variances are indexed by.
/// * `variance` - One vector of per-gene variance per group.
/// * `hvg` - The selected gene indices.
///
/// ### Returns
///
/// A list with `genes`, a genes by groups `variance` matrix and `hvg`.
pub(crate) fn residual_variance_to_r_list(
    genes: &[usize],
    variance: &[Vec<f64>],
    hvg: &[usize],
) -> List {
    let n_genes = genes.len();
    let n_groups = variance.len();

    // Column-major, which is what R expects of a matrix, and the groups are
    // the columns so a single-group fit comes back as one column.
    let mut flat: Vec<f64> = Vec::with_capacity(n_genes * n_groups);
    for group in variance {
        flat.extend_from_slice(group);
    }

    let variance_matrix =
        RArray::new_matrix(n_genes, n_groups, |row, col| flat[col * n_genes + row]);

    list!(
        genes = genes.iter().map(|&g| g as i32).collect::<Vec<i32>>(),
        variance = variance_matrix,
        hvg = hvg.iter().map(|&g| g as i32).collect::<Vec<i32>>()
    )
}

/////////
// PCA //
/////////

/// Calculates PCA on Pearson residuals for single cell
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Runs the PCA on the residuals a fitted model implies, rather than on the
/// stored normalised layer. The residual columns are dense by construction,
/// since a zero count still has a residual, so there is no sparse or streaming
/// variant of this path.
///
/// Two settings are refused rather than ignored: the `PFlogPF` transform, which
/// belongs to the normalised layer, and variance normalisation, which would
/// flatten the very ranking the residuals produce.
///
/// @param f_path_gene String. Path to the `counts_genes.bin` file.
/// @param residual_fit List. A fit from [bixverse::rs_sc_fit_residuals()].
/// @param no_pcs Integer. Number of PCs to calculate.
/// @param pca_params Named list. Contains the parameters to use for this PCA
/// run. `clr` and `normalise_variance` must both be `FALSE`.
/// @param cell_indices Integer vector. The cell indices to use. (0-indexed!)
/// Must be the selection the fit was fitted on.
/// @param gene_indices Integer vector. The gene indices to use. (0-indexed!)
/// Every one must be covered by the fit.
/// @param seed Integer. Random seed for the randomised SVD.
/// @param return_scaled Boolean. Shall the scaled data be returned.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @returns A list with the following items
/// \itemize{
///   \item scores - The samples projected on the PCA space.
///   \item loadings - The loadings of the features for the PCA.
///   \item singular_values - The singular values for the PCA.
///   \item scaled - The scaled matrix if `return_scaled = TRUE`, otherwise
///   `NULL`.
/// }
///
/// @export
///
/// @keywords internal
#[allow(clippy::too_many_arguments)]
#[extendr]
fn rs_sc_pca_residuals(
    f_path_gene: &str,
    residual_fit: List,
    no_pcs: usize,
    pca_params: List,
    cell_indices: Vec<i32>,
    gene_indices: Vec<i32>,
    seed: usize,
    return_scaled: bool,
    verbose: usize,
) -> Result<List, extendr_api::Error> {
    let cell_set = cell_indices.r_int_convert();
    let gene_set = gene_indices.r_int_convert();

    let pca_params = SingleCellPcaParams::from_r_list(pca_params)?;

    let reader = ParallelSparseReader::new(f_path_gene).to_extendr()?;
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
            return_scaled,
            verbose,
        )
    })?;

    let singular_values_f64: Vec<f64> = res.2.iter().map(|&x| x as f64).collect();
    let scaled = res.3.map(|s| faer_to_r_matrix(s.as_ref()));

    Ok(list!(
        scores = faer_to_r_matrix(res.0.as_ref()),
        loadings = faer_to_r_matrix(res.1.as_ref()),
        singular_values = singular_values_f64,
        scaled = scaled
    ))
}

//////////////////////
// Corrected counts //
//////////////////////

/// Writes scTransform-corrected counts to a new store
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Reverses the residual transform with every latent variable, the library size
/// included, held at its median, so the depth structure is removed while the
/// per-sample intercept is kept. The result is written as a new gene-major
/// store.
///
/// The output is re-indexed: its gene axis is the model's, so gene `j` in the
/// written store is `genes[j + 1]` of the source. It also has no normalised
/// layer, since corrected counts carry no library size to scale to.
///
/// @param f_path_gene String. Path to the `counts_genes.bin` file.
/// @param residual_fit List. A scTransform fit from
/// [bixverse::rs_sc_fit_residuals()]. The analytic Pearson model has no
/// corrected-count equivalent.
/// @param cell_indices Integer vector. The cell indices to use. (0-indexed!)
/// Must be the selection the fit was fitted on.
/// @param f_path_out String. Path of the gene-major file to write.
/// @param gene_batch_size Integer or `NULL`. Genes held in memory per batch.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @returns A list with the following items
/// \itemize{
///   \item f_path - The file that was written.
///   \item genes - The source gene indices of the written axis. (0-indexed!)
///   \item n_genes - Number of genes written.
///   \item n_cells - Number of cells written.
/// }
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_sct_corrected_counts(
    f_path_gene: &str,
    residual_fit: List,
    cell_indices: Vec<i32>,
    f_path_out: &str,
    gene_batch_size: Option<usize>,
    verbose: usize,
) -> Result<List, extendr_api::Error> {
    let cell_set = cell_indices.r_int_convert();
    let opts = stream_opts(gene_batch_size, SctStreamOpts::default().seed, verbose);

    let reader = ParallelSparseReader::new(f_path_gene).to_extendr()?;
    let fit = parse_residual_fit(residual_fit, &cell_set)?;

    let RResidualFit::Sct {
        fit,
        log10_umi,
        covariates,
    } = &fit
    else {
        return Err(Error::Other(
            "Corrected counts need a scTransform fit; the analytic Pearson model \
             has no corrected-count equivalent."
                .to_string(),
        ));
    };

    let cells = SctCellContext::new(log10_umi, covariates).to_extendr()?;
    let source = SctResiduals::new(&fit.models, cells, fit.group_of_cell.clone()).to_extendr()?;

    sct_corrected_counts(&reader, &source, &cell_set, f_path_out, opts).to_extendr()?;

    let genes = source.genes();

    Ok(list!(
        f_path = f_path_out,
        genes = genes.iter().map(|&g| g as i32).collect::<Vec<i32>>(),
        n_genes = genes.len() as i32,
        n_cells = cell_set.len() as i32
    ))
}

/// Rebuilds the cell-major companion of a gene-major store
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Writes the `counts_cells.bin` twin of a `counts_genes.bin` file. Memory is
/// bounded by phasing over cells: each phase holds one window of cells and
/// re-reads the gene file to fill it, so peak memory is the window rather than
/// the matrix.
///
/// @param f_path_in String. Path to the gene-major source file.
/// @param f_path_out String. Path of the cell-major file to write.
/// @param cells_per_phase Integer. Cells held in memory at once.
/// @param gene_batch_size Integer. Genes read per batch within a phase.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @returns String. The path that was written.
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_sc_gene_store_to_cell_store(
    f_path_in: &str,
    f_path_out: &str,
    cells_per_phase: usize,
    gene_batch_size: usize,
    verbose: usize,
) -> Result<String, extendr_api::Error> {
    gene_store_to_cell_store(
        f_path_in,
        f_path_out,
        cells_per_phase,
        gene_batch_size,
        verbose,
    )
    .to_extendr()?;

    Ok(f_path_out.to_string())
}
