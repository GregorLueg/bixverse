use faer::{Mat, MatRef};
use manifolds_rs::{umap, UmapParams};

/// Wrapper function into the UMAP implementation in `manifolds-rs`
///
/// ### Params
///
/// ### Returns
pub fn bx_umap(
    data: MatRef<f32>,
    k: usize,
    n_dim: usize,
    ann_type: String,
    optim: String,
    seed: usize,
    verbose: bool,
) -> Mat<f32> {
    let res = umap(
        data,
        n_dim,
        k,
        optim,
        ann_type,
        &UmapParams::default(),
        None,
        None,
        seed,
        verbose,
    );

    let ncol = res.len();
    let nrow = res[0].len();

    Mat::from_fn(nrow, ncol, |i, j| res[j][i])
}
