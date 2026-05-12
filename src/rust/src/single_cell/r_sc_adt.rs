use bixverse_rs::prelude::*;
use bixverse_rs::single_cell::multi_modal::dsb::*;
use extendr_api::Nullable::Null;
use extendr_api::*;
use rayon::prelude::*;

////////////////////
// extendr Module //
////////////////////

extendr_module! {
    mod r_sc_adt;
    // processing
    fn rs_adt_clr;
    fn rs_dsb;
}

////////////////
// Processing //
////////////////

/// Applies CLR normalisation on ADT counts (Seurat-style, per cell)
///
/// @param counts R matrix of shape cells x features.
///
/// @returns CLR-transformed matrix.
///
/// @export
#[extendr]
fn rs_adt_clr(counts: RMatrix<f64>) -> RMatrix<f64> {
    let counts = r_matrix_to_faer(&counts);
    let nrow = counts.nrows();
    let ncol = counts.ncols();

    // per-row geometric-mean factor
    let g: Vec<f64> = (0..nrow)
        .into_par_iter()
        .map(|i| {
            let mut s = 0.0_f64;
            for j in 0..ncol {
                let v = counts[(i, j)];
                if v > 0.0 {
                    s += v.ln_1p();
                }
            }
            (s / ncol as f64).exp()
        })
        .collect();

    RMatrix::new_matrix(nrow, ncol, |i, j| (counts[(i, j)] / g[i]).ln_1p())
}

/// Run DSB normalisation on raw ADT counts
///
/// @description This function applies the DSB algorithm to normalise
/// antibody-derived tag (ADT) counts from CITE-seq experiments. Two variants
/// are supported. When `background_counts` is provided, per-protein ambient
/// background is estimated from empty droplets ("Step I" of the original
/// paper). When `background_counts` is `NULL`, per-protein background is
/// estimated by a two-component k-means on the log-transformed cell counts,
/// with the lower centroid taken as the background level. An optional second
/// step removes cell-to-cell technical noise by regressing out PC1 of a noise
/// matrix built from isotype controls (if available) and the per-cell
/// background mean.
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
/// @param seed Integer. Random seed for k-means initialisation.
/// @param verbose Boolean. Print progress messages.
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
///   technical component regressed out in Step II.
///   \item cellwise_background_mean - Numeric vector of length `n_cells`,
///   or empty vector if `dsb_params$denoise_counts = FALSE`. Per-cell
///   background mean from the 2-component k-means clustering.
/// }
///
/// @export
#[extendr]
fn rs_dsb(
    raw_counts: RMatrix<f64>,
    background_counts: Nullable<RMatrix<f64>>,
    isotype_indices: Vec<i32>,
    dsb_params: List,
    scale_factor: String,
    seed: usize,
    verbose: bool,
) -> Result<List> {
    let raw_counts = r_matrix_to_faer_fp32(&raw_counts);
    let dsb_params = DsbParams::from_r_list(dsb_params, isotype_indices)?;

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
