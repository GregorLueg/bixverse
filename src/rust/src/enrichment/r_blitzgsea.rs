use bixverse_rs::enrichment::blitzgsea::{BlitzGseaParams, blitzgsea_score, calibrate_null};
use bixverse_rs::enrichment::enrichment_r_wrapper::{
    blitzgsea_null_from_list, blitzgsea_null_to_list,
};
use bixverse_rs::prelude::*;
use extendr_api::prelude::*;

/////////////
// extendR //
/////////////

extendr_module! {
    mod r_blitzgsea;
    fn rs_blitzgsea_calibrate;
    fn rs_blitzgsea_score;
}

//////////////////////
// Helper functions //
//////////////////////

/// Pull the gene set index vectors out of an R list
///
/// ### Params
///
/// * `pathways` - R list holding one integer vector of index positions per gene
///   set, as [`rs_get_gs_indices`](crate::enrichment::r_gsea) returns them.
///
/// ### Returns
///
/// The index vectors, or an error naming the offending gene set.
fn r_list_to_index_vec(pathways: List) -> extendr_api::Result<Vec<Vec<i32>>> {
    (0..pathways.len())
        .map(|i| {
            pathways.elt(i)?.as_integer_vector().ok_or_else(|| {
                extendr_api::Error::Other(format!(
                    "Gene set {} is not an integer vector of index positions.",
                    i + 1
                ))
            })
        })
        .collect()
}

////////////////////
// Core functions //
////////////////////

/// Calibrate the blitzGSEA gamma null for a signature
///
/// @description
/// `r lifecycle::badge("experimental")`
///
/// Draws random gene sets across a log-spaced grid of anchor sizes, fits gamma
/// tails to the resulting enrichment scores and smooths the fitted parameters
/// across sizes. Nothing about any gene set library enters here, so one
/// calibration serves every library scored against that signature.
///
/// @param stats Numeric vector. The gene level statistic. Needs to be sorted in
/// descending nature.
/// @param blitz_params List. The blitzGSEA parameters, see
/// [bixverse::params_blitzgsea()]. Recognised elements are `permutations`,
/// `anchors`, `symmetric`, `centre`, `ks_test` and `seed`; anything else is
/// ignored and any missing element takes its default.
///
/// @return List with the following elements
/// \itemize{
///     \item anchor_sizes Numeric vector. The anchor set sizes, ascending.
///     \item shape_pos Numeric vector. Smoothed positive-tail gamma shape.
///     \item scale_pos Numeric vector. Smoothed positive-tail gamma scale.
///     \item shape_neg Numeric vector. Smoothed negative-tail gamma shape.
///     \item scale_neg Numeric vector. Smoothed negative-tail gamma scale.
///     \item pos_ratio Numeric vector. Smoothed fraction of positive null
///     scores at each anchor.
///     \item ks_pos Float. Mean goodness-of-fit p-value for the positive tail.
///     \item ks_neg Float. Mean goodness-of-fit p-value for the negative tail.
///     \item centred Boolean. Whether the signature was centred.
/// }
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_blitzgsea_calibrate(stats: &[f64], blitz_params: List) -> extendr_api::Result<List> {
    let params = BlitzGseaParams::from_r_list(blitz_params)?;

    let null = calibrate_null(stats, &params).to_extendr()?;

    Ok(blitzgsea_null_to_list(&null))
}

/// Score gene sets against a calibrated blitzGSEA null
///
/// @description
/// `r lifecycle::badge("experimental")`
///
/// Each gene set costs one enrichment score plus one gamma tail evaluation.
/// The gene sets are expected to have been filtered to the desired size bounds
/// and intersected with the signature already.
///
/// @param stats Numeric vector. The gene level statistic. Needs to be sorted in
/// descending nature and be the same signature the null was calibrated on.
/// @param pathways List. One integer vector of index positions per gene set,
/// indexed to R's 1-indexing. Order and duplicates do not matter.
/// @param null_model List. The calibrated null from
/// [bixverse::rs_blitzgsea_calibrate()].
/// @param blitz_params List. The blitzGSEA parameters, see
/// [bixverse::params_blitzgsea()]. Only `centre` is read here and it has to
/// match what the calibration used.
///
/// @return List with the following elements
/// \itemize{
///     \item es Numeric vector. Enrichment scores for the gene sets.
///     \item nes Numeric vector. Normalised enrichment scores.
///     \item pvals Numeric vector. Two-sided p-values from the gamma
///     approximation.
///     \item sidak Numeric vector. Sidak-adjusted p-values.
///     \item fdr Numeric vector. Benjamini-Hochberg adjusted p-values.
///     \item size Integer vector. Gene set size after intersection.
///     \item leading_edge List of integer vectors with the leading edge index
///     positions, indexed to R's 1-indexing.
/// }
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_blitzgsea_score(
    stats: &[f64],
    pathways: List,
    null_model: List,
    blitz_params: List,
) -> extendr_api::Result<List> {
    let params = BlitzGseaParams::from_r_list(blitz_params)?;
    let null = blitzgsea_null_from_list(null_model)?;
    let pathways = r_list_to_index_vec(pathways)?;

    let res = blitzgsea_score(stats, &pathways, &null, &params, true).to_extendr()?;

    let mut leading_edge = List::new(res.leading_edge.len());
    for (i, x) in res.leading_edge.iter().enumerate() {
        leading_edge.set_elt(i, Robj::from(x.clone()))?;
    }

    Ok(list!(
        es = res.es,
        nes = res.nes,
        pvals = res.pvals,
        sidak = res.sidak,
        fdr = res.fdr,
        size = res.size,
        leading_edge = leading_edge
    ))
}
