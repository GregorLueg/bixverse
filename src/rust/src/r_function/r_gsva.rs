use extendr_api::prelude::*;

use crate::helpers::gsva::{get_gsva_gs_indices, gsva};
use crate::utils::r_rust_interface::{faer_to_r_matrix, NamedMatrix};

/// @export
#[extendr]
fn rs_gsva(
    exp: RMatrix<f64>,
    gs_list: List,
    tau: f64,
    gaussian: bool,
    max_diff: bool,
    abs_rank: bool,
) -> extendr_api::Result<extendr_api::RArray<f64, [usize; 2]>> {
    let exp = NamedMatrix::new(&exp);

    let gs_indices = get_gsva_gs_indices(&exp, gs_list)?;

    let results = gsva(&exp.values, &gs_indices, gaussian, tau, max_diff, abs_rank);

    Ok(faer_to_r_matrix(results.as_ref()))
}

extendr_module! {
    mod r_gsva;
    fn rs_gsva;
}
