use extendr_api::*;

use crate::single_cell::methods::umap::*;
use crate::utils::r_rust_interface::*;

extendr_module! {
    mod r_sc_embd;
    fn rs_umap;
}

/// @export
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_umap(
    embd: RMatrix<f64>,
    n_dim: usize,
    min_dist: f64,
    spread: f64,
    ann_type: String,
    optim: String,
    init: String,
    k: usize,
    seed: usize,
    verbose: bool,
) -> RMatrix<f64> {
    let embd = r_matrix_to_faer_fp32(&embd);

    let res = bxv_umap(
        embd.as_ref(),
        k,
        min_dist as f32,
        spread as f32,
        n_dim,
        ann_type,
        optim,
        init,
        seed,
        verbose,
    );

    faer_to_r_matrix(res.as_ref())
}
