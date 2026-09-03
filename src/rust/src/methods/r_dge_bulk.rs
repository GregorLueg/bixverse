//! Bindings for edgeR's quasi-likelihood chain, see Chen, Lun and Smyth,
//! F1000Research, 2016.
//!
//! The numerics are `edge-rs`, assembled into the edgeR chain in
//! [`bixverse_rs::methods::dge_bulk`]. Anything with a counts matrix of the
//! tested axis by samples goes through here, so both plain bulk and Milo's
//! neighbourhood counts use this one binding.

use bixverse_rs::methods::dge_bulk::{run_edger_ql, EdgeRQlParams};
use bixverse_rs::methods::methods_r_wrapper::TestedFromR;
use bixverse_rs::prelude::*;
use edge_rs::glm::test::Tested;
use extendr_api::*;
use std::collections::HashMap;

/////////////
// extendR //
/////////////

extendr_module! {
    mod r_dge_bulk;
    fn rs_edger_ql;
}

///////////////
// Functions //
///////////////

/// Run the edgeR quasi-likelihood chain on a count matrix
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Runs `filterByExpr` -> `calcNormFactors` -> `glmQLFit` -> `glmQLFTest`,
/// implemented in Rust via the `edge-rs` crate and gated against edgeR 4.8.2.
/// The tested axis does not have to be genes: Milo's neighbourhood counts are
/// tested with the same call, with `filter = FALSE`.
///
/// @param counts Numeric matrix. Raw counts of features x samples. Must not
/// be normalised or log-transformed.
/// @param design Numeric matrix. The design matrix of samples x coefficients,
/// including the intercept. Needs at least two columns, since the null model
/// has to retain one.
/// @param edger_params Named list. The edgeR parameters, see
/// [bixverse::params_edger_ql()], plus either `coef` (0-indexed(!) design
/// columns to drop from the null model) or `contrast` (column-major weights
/// with `n_contrasts` columns).
///
/// @return A list with the following elements
/// \itemize{
///   \item features_to_keep - Boolean. Which features survived the filters.
///   Spans the full feature axis of `counts`.
///   \item log_fc - Log2 fold changes of the tested coefficient or contrast.
///   \item log_cpm - Average log2 counts per million.
///   \item f_stat - The quasi-likelihood F statistic.
///   \item p_values - Raw p-values.
///   \item fdr - Benjamini-Hochberg adjusted p-values.
/// }
///
/// @references Chen, Lun and Smyth, F1000Research, 2016
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_edger_ql(counts: RMatrix<f64>, design: RMatrix<f64>, edger_params: List) -> Result<List> {
    let n_features = counts.nrows();
    let n_samples = counts.ncols();
    let n_coef = design.ncols();

    if design.nrows() != n_samples {
        return Err(Error::Other(format!(
            "The design has {} rows against {} samples in the counts.",
            design.nrows(),
            n_samples
        )));
    }

    let params_map: HashMap<&str, Robj> = r_list_to_map(edger_params.clone())?;
    let tested = Tested::from_r_map(&params_map)?;
    let params = EdgeRQlParams::from_r_list(edger_params)?;

    // Both are column-major coming out of R and `run_edger_ql` reads them
    // row-major.
    let counts = mat_to_flat_row_major(r_matrix_to_faer(&counts));
    let design = mat_to_flat_row_major(r_matrix_to_faer(&design));

    let res = run_edger_ql(
        &counts, n_features, n_samples, &design, n_coef, &tested, &params,
    )
    .to_extendr()?;

    Ok(list!(
        features_to_keep = res.genes_to_keep,
        log_fc = res.log_fc,
        log_cpm = res.log_cpm,
        f_stat = res.f_stat,
        p_values = res.p_val,
        fdr = res.fdr
    ))
}
