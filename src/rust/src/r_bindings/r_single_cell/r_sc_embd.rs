use extendr_api::*;

use crate::single_cell::methods::umap::*;
use crate::utils::r_rust_interface::*;

extendr_module! {
    mod r_sc_embd;
    fn rs_umap;
}

/// @export
#[extendr]
fn rs_umap(
    embd: RMatrix<f64>,
    n_dim: usize,
    ann_type: String,
    optim: String,
    k: usize,
    seed: usize,
    verbose: bool,
) -> RMatrix<f64> {
    let embd = r_matrix_to_faer_fp32(&embd);

    let res = bx_umap(embd.as_ref(), k, n_dim, ann_type, optim, seed, verbose);

    faer_to_r_matrix(res.as_ref())
}
