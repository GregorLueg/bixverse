use bixverse_r::enrichment::lm_activity::{lm_activity, LmActivityParams};
use bixverse_r::enrichment_r_wrapper::r_list_to_feature_set_weights;
use bixverse_r::prelude::*;
use extendr_api::prelude::*;

/////////////
// extendR //
/////////////

extendr_module! {
    mod r_lm_activity;
    fn rs_lm_activity;
}

////////////////////
// Core functions //
////////////////////

/// ULM / MLM activity scores (decoupler-style)
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Fits a univariate (ULM) or multivariate (MLM) linear model per observation
/// of the feature weights of each set against the observed values, and returns
/// the t-value or the coefficient as the activity score. Observations are the
/// columns of `mat`.
///
/// @param mat Numerical matrix. Features x observations. Must not contain
/// non-finite values.
/// @param gs_list List. One element per set, each a list with `idx` (integer,
/// 0-indexed feature positions) and optionally `weight` (numeric, aligned with
/// `idx`; absent means weight 1).
/// @param lm_params List. With `method` (`"ulm"` or `"mlm"`), `stat` (`"t"` or
/// `"coef"`) and `p_values` (boolean).
///
/// @returns A list with the following elements
/// \itemize{
///   \item es - Numerical matrix of observations x sets with the scores.
///   \item p_values - Numerical matrix of observations x sets with the
///   two-sided p-values, or `NULL` if not requested.
/// }
///
/// @export
///
/// @keywords internal
///
/// @references Badia-i-Mompel, et al., Bioinform. Adv., 2022.
#[extendr]
fn rs_lm_activity(
    mat: RMatrix<f64>,
    gs_list: List,
    lm_params: List,
) -> extendr_api::Result<List> {
    let mat = r_matrix_to_faer(&mat);
    let net = r_list_to_feature_set_weights(gs_list, mat.nrows())?;
    let params = LmActivityParams::from_r_list(lm_params)?;

    let res = lm_activity(mat, &net, Some(params)).to_extendr()?;

    let p_values = match res.p_values {
        Some(p) => faer_to_r_matrix(p.as_ref()).into_robj(),
        None => NULL.into_robj(),
    };

    Ok(list!(es = faer_to_r_matrix(res.es.as_ref()), p_values = p_values))
}
