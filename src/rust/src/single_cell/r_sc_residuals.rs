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

use bixverse_r::prelude::*;
use bixverse_r::single_cell::sc_data::bin_merge_io::gene_store_to_cell_store;
use bixverse_r::single_cell::sc_processing::analytic_pearson::model::AprParams;
use bixverse_r::single_cell::sc_processing::pca::SingleCellPcaParams;
use bixverse_r::single_cell::sc_processing::residual_fit::{
    fit_residuals_on_store, pca_residuals_on_store, residual_variance_on_store,
    sct_corrected_counts_on_store, ResidualFit, ResidualModel, ResidualParams,
    ResidualVarianceRes,
};
use bixverse_r::single_cell::sc_processing::sctransform::model::{SctCovariates, SctParams};
use bixverse_r::single_cell::sc_processing::sctransform::stream::{SctGroupedFit, SctStreamOpts};
use bixverse_r::single_cell::sc_processing::analytic_pearson::stream::AprGroupedFit;
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

//////////////
// Fit <> R //
//////////////

/// Serialise a fitted residual model for R.
///
/// ### Params
///
/// * `fit` - The fitted model and the cells it was fitted on.
///
/// ### Returns
///
/// The list [`parse_residual_fit`] reads back. The fitted cells are echoed so
/// a later call can check it is being handed the same selection; the
/// analytic Pearson list also carries the shared gene axis, which R wants for
/// the HVG bookkeeping.
pub(crate) fn residual_fit_to_r_list(fit: &ResidualFit) -> List {
    let echoed: Vec<i32> = fit.cell_indices.iter().map(|&c| c as i32).collect();

    match &fit.model {
        ResidualModel::Sct {
            fit: grouped,
            log10_umi,
            covariates,
        } => extend_list(
            grouped.to_r_list(),
            &[
                ("method", Robj::from(METHOD_SCT)),
                ("cell_indices", Robj::from(echoed)),
                ("covariates", Robj::from(covariates.to_r_list())),
                ("log10_umi", Robj::from(log10_umi.clone())),
            ],
        ),
        ResidualModel::Apr { fit: grouped } => {
            let genes: Vec<i32> = fit.genes().iter().map(|&g| g as i32).collect();
            extend_list(
                grouped.to_r_list(),
                &[
                    ("method", Robj::from(METHOD_APR)),
                    ("cell_indices", Robj::from(echoed)),
                    ("genes", Robj::from(genes)),
                ],
            )
        }
    }
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
/// Does not check the cells: every crate entry point taking a [`ResidualFit`]
/// refuses a selection other than the fitted one.
///
/// ### Params
///
/// * `fit` - The serialised fit.
///
/// ### Returns
///
/// The [`ResidualFit`], or an error naming what is missing.
pub(crate) fn parse_residual_fit(fit: List) -> Result<ResidualFit, extendr_api::Error> {
    let method: String = fit
        .index("method")
        .map_err(|_| Error::Other("The residual fit is missing 'method'".to_string()))?
        .as_str()
        .ok_or_else(|| Error::Other("The residual fit 'method' is not a string".to_string()))?
        .to_string();

    let cell_indices: Vec<usize> = fit
        .index("cell_indices")
        .map_err(|_| Error::Other("The residual fit is missing 'cell_indices'".to_string()))?
        .as_integer_slice()
        .ok_or_else(|| Error::Other("The residual fit 'cell_indices' is not integer".to_string()))?
        .iter()
        .map(|&c| c as usize)
        .collect();

    let model = match method.as_str() {
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

            ResidualModel::Sct {
                fit: SctGroupedFit::from_r_list(fit)?,
                log10_umi,
                covariates,
            }
        }
        METHOD_APR => ResidualModel::Apr {
            fit: AprGroupedFit::from_r_list(fit)?,
        },
        other => {
            return Err(Error::Other(format!(
                "Unknown residual method '{other}'. Expected '{METHOD_SCT}' or '{METHOD_APR}'."
            )))
        }
    };

    Ok(ResidualFit {
        model,
        cell_indices,
    })
}

/// The residual method and its parameters from the R arguments.
///
/// ### Params
///
/// * `method` - `"sctransform"` or `"analytic_pearson"`.
/// * `params` - The method's parameters as an R list.
///
/// ### Returns
///
/// The [`ResidualParams`], or an error on an unknown method.
pub(crate) fn residual_params(
    method: &str,
    params: List,
) -> Result<ResidualParams, extendr_api::Error> {
    match method {
        METHOD_SCT => Ok(ResidualParams::Sct(SctParams::from_r_list(params)?)),
        METHOD_APR => Ok(ResidualParams::Apr(AprParams::from_r_list(params)?)),
        other => Err(Error::Other(format!(
            "Unknown residual method '{other}'. Expected '{METHOD_SCT}' or '{METHOD_APR}'."
        ))),
    }
}

/// The covariates of a fit; only scTransform reads them.
///
/// ### Params
///
/// * `params` - The resolved method parameters.
/// * `covariates` - The covariates as an R list.
///
/// ### Returns
///
/// The parsed covariates for scTransform, empty ones otherwise.
pub(crate) fn residual_covariates(
    params: &ResidualParams,
    covariates: List,
) -> Result<SctCovariates, extendr_api::Error> {
    match params {
        ResidualParams::Sct(_) => SctCovariates::from_r_list(covariates),
        ResidualParams::Apr(_) => Ok(SctCovariates::default()),
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
pub(crate) fn stream_opts(
    gene_batch_size: Option<usize>,
    seed: u64,
    verbose: usize,
) -> SctStreamOpts {
    let defaults = SctStreamOpts::default();

    SctStreamOpts {
        gene_batch_size: gene_batch_size.or(defaults.gene_batch_size),
        seed,
        verbose,
    }
}

/// Group labels for a fit.
///
/// ### Params
///
/// * `group_of_cell` - The 0-based labels from R, or `None` for one model.
/// * `n_cells` - Number of selected cells.
///
/// ### Returns
///
/// One label per selected cell (or `None`), or an error when the length
/// disagrees or a label is negative.
pub(crate) fn resolve_groups(
    group_of_cell: Option<Vec<i32>>,
    n_cells: usize,
) -> Result<Option<Vec<u32>>, extendr_api::Error> {
    let Some(groups) = group_of_cell else {
        return Ok(None);
    };
    if groups.len() != n_cells {
        return Err(Error::Other(format!(
            "'group_of_cell' has {} entries but {n_cells} cells were selected",
            groups.len()
        )));
    }
    groups
        .into_iter()
        .map(|g| {
            u32::try_from(g)
                .map_err(|_| Error::Other("'group_of_cell' has a negative entry".to_string()))
        })
        .collect::<Result<Vec<u32>, _>>()
        .map(Some)
}

/// Pack the residual variance sweep into the list R expects.
///
/// ### Params
///
/// * `res` - The per-group variances and the selected genes.
///
/// ### Returns
///
/// A list with `genes`, a genes by groups `variance` matrix and `hvg`.
pub(crate) fn residual_variance_to_r_list(res: &ResidualVarianceRes) -> List {
    let n_genes = res.genes.len();
    let n_groups = res.variance.len();

    // groups are the columns, so a single-group fit comes back as one column
    let variance_matrix =
        RArray::new_matrix(n_genes, n_groups, |row, col| res.variance[col][row]);

    list!(
        genes = res.genes.iter().map(|&g| g as i32).collect::<Vec<i32>>(),
        variance = variance_matrix,
        hvg = res.hvg.iter().map(|&g| g as i32).collect::<Vec<i32>>()
    )
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
/// library sizes, the per-cell totals of the analytic Pearson fit and, for small
/// cell selections, the counts themselves instead of sweeping the gene file.
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
    let groups = resolve_groups(group_of_cell, cell_set.len())?;
    let params = residual_params(method, params)?;
    let covariates = residual_covariates(&params, covariates)?;

    let fit = fit_residuals_on_store(
        f_path_gene,
        f_path_cell,
        &cell_set,
        groups.as_deref(),
        covariates,
        params,
        stream_opts(gene_batch_size, seed, verbose),
    )
    .to_extendr()?;

    Ok(residual_fit_to_r_list(&fit))
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
/// @param f_path_cell String. Path to the `counts_cells.bin` file. Small cell
/// selections are read from here instead of sweeping the gene file.
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
    f_path_cell: &str,
    residual_fit: List,
    cell_indices: Vec<i32>,
    n_hvg: usize,
    gene_batch_size: Option<usize>,
    verbose: usize,
) -> Result<List, extendr_api::Error> {
    let fit = parse_residual_fit(residual_fit)?;

    let res = residual_variance_on_store(
        f_path_gene,
        f_path_cell,
        &fit,
        &cell_indices.r_int_convert(),
        n_hvg,
        stream_opts(gene_batch_size, SctStreamOpts::default().seed, verbose),
    )
    .to_extendr()?;

    Ok(residual_variance_to_r_list(&res))
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
/// @param f_path_cell String. Path to the `counts_cells.bin` file. Small cell
/// selections are read from here instead of sweeping the gene file.
/// @param residual_fit List. A fit from [bixverse::rs_sc_fit_residuals()].
/// @param no_pcs Integer. Number of PCs to calculate.
/// @param pca_params Named list. Contains the parameters to use for this PCA
/// run. `clr` and `normalise_variance` must both be `FALSE`.
/// @param cell_indices Integer vector. The cell indices to use. (0-indexed!)
/// Must be the selection the fit was fitted on.
/// @param gene_indices Integer vector. The gene indices to use. (0-indexed!)
/// Every one must be covered by the fit.
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
fn rs_sc_pca_residuals(
    f_path_gene: &str,
    f_path_cell: &str,
    residual_fit: List,
    no_pcs: usize,
    pca_params: List,
    cell_indices: Vec<i32>,
    gene_indices: Vec<i32>,
    seed: usize,
    verbose: usize,
) -> Result<List, extendr_api::Error> {
    let pca_params = SingleCellPcaParams::from_r_list(pca_params)?;
    let fit = parse_residual_fit(residual_fit)?;

    let (scores, loadings, singular_values) = pca_residuals_on_store(
        f_path_gene,
        f_path_cell,
        &fit,
        &cell_indices.r_int_convert(),
        &gene_indices.r_int_convert(),
        no_pcs,
        &pca_params,
        seed,
        verbose,
    )
    .to_extendr()?;

    Ok(list!(
        scores = faer_to_r_matrix(scores.as_ref()),
        loadings = faer_to_r_matrix(loadings.as_ref()),
        singular_values = singular_values.r_float_convert()
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
    let fit = parse_residual_fit(residual_fit)?;

    let genes = sct_corrected_counts_on_store(
        f_path_gene,
        &fit,
        &cell_set,
        f_path_out,
        stream_opts(gene_batch_size, SctStreamOpts::default().seed, verbose),
    )
    .to_extendr()?;

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
