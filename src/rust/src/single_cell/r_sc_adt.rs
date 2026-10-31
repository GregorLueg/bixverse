use bixverse_r::mm_r_wrappers::dsb_params_from_r_list;
use bixverse_r::prelude::*;
use bixverse_r::single_cell::multi_modal::adt::clr::adt_clr;
use bixverse_r::single_cell::multi_modal::adt::dsb::*;
use bixverse_r::single_cell::sc_data::h5_10x_adt_io::{read_tenx_h5_modality, TenxDenseModality};
use bixverse_r::single_cell::sc_data::h5_10x_io::parse_tenx_version;
use extendr_api::Nullable::Null;
use extendr_api::*;


////////////////////
// extendr Module //
////////////////////

extendr_module! {
    mod r_sc_adt;
    // i/o
    fn rs_read_tenx_h5_modality;
    // processing
    fn rs_adt_clr;
    fn rs_dsb;
}

/////////
// I/O //
/////////

/// Loads in a modality from a 10x h5 file
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Reads one feature type (e.g. `"Antibody Capture"`) from a 10x h5 file as a
/// dense matrix. Needs a v3 file; v2 files carry no feature types and error.
///
/// @param f_path String. The path to the h5 file.
/// @param version String. The 10x version. If `"auto"`, the version is
/// detected from the file.
/// @param feature_type String. The feature type to return. Matched against
/// the whitespace-trimmed feature types of the file.
///
/// @returns A list with:
/// \itemize{
///   \item counts - Dense numerical matrix of cells x features.
///   \item barcodes - Character vector of the cell barcodes, in file order.
///   \item features - Character vector of the feature names.
/// }
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_read_tenx_h5_modality(f_path: String, version: String, feature_type: String) -> Result<List> {
    let tenx_version = match version.to_lowercase().as_str() {
        "auto" => None,
        other => parse_tenx_version(other),
    };

    let modality_res: TenxDenseModality =
        read_tenx_h5_modality(f_path, tenx_version, &feature_type).to_extendr()?;

    let counts = RMatrix::new_matrix(modality_res.n_cells, modality_res.n_features, |i, j| {
        modality_res.counts[j * modality_res.n_cells + i]
    });

    Ok(list!(
        counts = counts,
        barcodes = modality_res.barcodes,
        features = modality_res.features
    ))
}

////////////////
// Processing //
////////////////

/////////////
// Helpers //
/////////////

/// Applies CLR normalisation on ADT counts
///
/// @description
/// `r lifecycle::badge("experimental")`
///
/// Normalises each cell (row) separately.
///
/// @param counts Numerical matrix of shape cells x features.
/// @param seurat_clr Boolean. If `TRUE` uses the Seurat variant
/// `log1p(x / g)` (non-negative); if `FALSE` uses the proper CLR
/// `log1p(x) - mean(log1p(x))` (mean-centred, can be negative).
///
/// @returns Numerical matrix of cells x features with the CLR-transformed
/// values.
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_adt_clr(counts: RMatrix<f64>, seurat_clr: bool) -> RMatrix<f64> {
    let counts = r_matrix_to_faer(&counts);
    faer_to_r_matrix(adt_clr(counts, seurat_clr).as_ref())
}

/// Run DSB normalisation on raw ADT counts
///
/// @description
/// `r lifecycle::badge("experimental")`
/// This function applies the DSB algorithm to normalise antibody-derived tag
/// (ADT) counts from CITE-seq experiments. Two variants are supported. When
/// `background_counts` is provided, per-protein ambient background is estimated
/// from empty droplets ("Step I" of the original paper). When
/// `background_counts` is `NULL`, per-protein background is estimated by a
/// two-component k-means on the `log(x + pseudocount)` cell counts, with the
/// lower centroid taken as the background level. An optional second step
/// removes cell-to-cell technical noise by regressing out PC1 of a noise
/// matrix built from isotype controls (if used) and the per-cell background
/// mean; without isotype controls the per-cell background mean itself is
/// regressed out. Optional per-protein quantile clipping runs last.
///
/// @param raw_counts Numeric matrix. Cells x proteins matrix of raw ADT
/// counts.
/// @param background_counts Optional numeric matrix. Cells x proteins matrix
/// of empty-droplet ADT counts used to estimate per-protein ambient
/// background. Must have the same number of columns as `raw_counts`. If
/// `NULL`, the function falls back to k-means-based background estimation.
/// @param isotype_indices Integer vector. 0-based column indices into
/// `raw_counts` identifying isotype control proteins. Used in Step II if
/// `dsb_params$use_isotype_controls = TRUE`. Pass an empty integer vector
/// if there are no isotype controls.
/// @param dsb_params List. DSB parameters.
/// @param scale_factor String. One of `"standardise"` or `"mean_subtract"`.
/// Only used when `background_counts` is provided. `"standardise"` subtracts
/// the per-protein background mean and divides by the per-protein background
/// SD. `"mean_subtract"` subtracts the mean only.
/// @param seed Integer. Random seed for the k-means initialisations.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @returns A list with the following elements:
/// \itemize{
///   \item norm_counts - Numeric matrix. The DSB-normalised cells x proteins
///   matrix.
///   \item protein_background_mean - Numeric vector of length `n_proteins`.
///   Per-protein background mean used in Step I.
///   \item protein_background_sd - Numeric vector of length `n_proteins`, or
///   empty vector if `background_counts` was `NULL`. Per-protein background
///   SD used in Step I.
///   \item technical_component - Numeric vector of length `n_cells`, or
///   empty vector if `dsb_params$denoise_counts = FALSE`. Per-cell
///   covariate regressed out in Step II: PC1 of the noise matrix, or the
///   per-cell background mean if no isotype controls are used.
///   \item cellwise_background_mean - Numeric vector of length `n_cells`,
///   or empty vector if `dsb_params$denoise_counts = FALSE`. Per-cell
///   background mean from a 2-component k-means on each cell's protein
///   vector.
/// }
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_dsb(
    raw_counts: RMatrix<f64>,
    background_counts: Nullable<RMatrix<f64>>,
    isotype_indices: Vec<i32>,
    dsb_params: List,
    scale_factor: String,
    seed: usize,
    verbose: usize,
) -> Result<List> {
    // prepare_background_source falls back to "standardise" on anything else
    if parse_scale_factor(&scale_factor).is_none() {
        return Err(Error::Other(format!(
            "Unknown scale_factor '{scale_factor}'. Use 'standardise' or 'mean_subtract'."
        )));
    }
    let raw_counts = r_matrix_to_faer_fp32(&raw_counts);
    let dsb_params = dsb_params_from_r_list(dsb_params, isotype_indices)?;

    let background_counts_provided = background_counts != Null;

    let background_counts = if background_counts_provided {
        let background_counts: RArray<f64, 2> = background_counts
            .into_robj()
            .as_matrix()
            .ok_or_else(|| Error::Other("'background_counts' is not a matrix".into()))?;

        let background_counts = r_matrix_to_faer_fp32(&background_counts);

        Some(background_counts)
    } else {
        None
    };

    let background_source =
        prepare_background_source(background_counts.as_ref().map(|m| m.as_ref()), scale_factor);

    let dsb_res: DsbResult = dsb_normalise(
        raw_counts.as_ref(),
        background_source,
        &dsb_params,
        verbose,
        seed,
    )
    .to_extendr()?;

    Ok(list!(
        norm_counts = faer_to_r_matrix(dsb_res.normalised.as_ref()),
        protein_background_mean = dsb_res.protein_background_mean,
        protein_background_sd = dsb_res.protein_background_sd.unwrap_or_default(),
        technical_component = dsb_res.technical_component.unwrap_or_default(),
        cellwise_background_mean = dsb_res.cellwise_background_mean.unwrap_or_default()
    ))
}
