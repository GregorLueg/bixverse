use extendr_api::prelude::*;

use crate::helpers::gsva::{get_gsva_gs_indices, gsva};
use crate::utils::r_rust_interface::{faer_to_r_matrix, NamedMatrix};

/// Rust version of the GSVA algorithm
///
/// @description
/// Rust-based implementation of the popular GSVA algorithm. Has further
/// performance optimisations compared to the original implementation.
///
/// @param exp Numerical matrix. The expression matrix with rows = genes, and
/// columns = samples
/// @param gs_list List. A list containing the pathway genes.
/// @param tau Tau parameter. Usual recommendation is to use `1.0` here. Larger
/// values emphasise the tails more.
/// @param gaussian If `TRUE` the Gaussian kernel will be used, if `FALSE` the
/// Poisson kernel will be used.
/// @param max_diff Scoring mode: `TRUE` = difference, `FALSE` = larger absolute
/// value
/// @param abs_rank If `TRUE` = pos-neg, `FALSE` = pos+neg
///
/// @return Returns a matrix of gene set ES scores x samples.
///
/// @export
#[extendr]
fn rs_gsva(
    exp: RMatrix<f64>,
    gs_list: List,
    tau: f64,
    gaussian: bool,
    max_diff: bool,
    abs_rank: bool,
    timings: bool,
) -> extendr_api::Result<extendr_api::RArray<f64, [usize; 2]>> {
    let exp = NamedMatrix::new(&exp);

    let gs_indices = get_gsva_gs_indices(&exp, gs_list)?;

    let results = gsva(
        &exp.values,
        &gs_indices,
        gaussian,
        tau,
        max_diff,
        abs_rank,
        timings,
    );

    Ok(faer_to_r_matrix(results.as_ref()))
}

extendr_module! {
    mod r_gsva;
    fn rs_gsva;
}
