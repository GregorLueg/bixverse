use extendr_api::List;
use faer::{Mat, MatRef};
use manifolds_rs::nearest_neighbours::NearestNeighbourParams;
use manifolds_rs::optimiser::OptimParams;
use manifolds_rs::{umap, UmapParams};

////////////
// Params //
////////////

/// Helper function to generate the NearestNeighbourParams for UMAP
///
/// This one is **specifically** designed for UMAP. Do not use anywere else
///
/// ### Params
///
/// * `nn_param_list` - The R list with the UMAP nearest neighbour parameters.
///
/// ### Returns
///
/// The `NearestNeighbourParams` for UMAP.
pub fn umap_nn_params(nn_param_list: List) -> NearestNeighbourParams<f32> {
    let nn_param_list = nn_param_list.into_hashmap();

    // general
    let dist_metric = nn_param_list
        .get("ann_dist")
        .and_then(|v| v.as_str())
        .unwrap_or("cosine")
        .to_string();

    // annoy
    let n_trees = nn_param_list
        .get("n_trees")
        .and_then(|v| v.as_integer())
        .unwrap_or(50) as usize;

    let search_budget = nn_param_list
        .get("search_budget")
        .and_then(|v| v.as_integer())
        .unwrap_or(2) as usize;

    // hnsw
    let m = nn_param_list
        .get("m")
        .and_then(|v| v.as_integer())
        .unwrap_or(16) as usize;

    let ef_construction = nn_param_list
        .get("ef_construction")
        .and_then(|v| v.as_integer())
        .unwrap_or(200) as usize;

    let ef_search = nn_param_list
        .get("ef_search")
        .and_then(|v| v.as_integer())
        .unwrap_or(2) as usize;

    // nndescent
    let max_iter = nn_param_list
        .get("max_iter")
        .and_then(|v| v.as_integer())
        .unwrap_or(25) as usize;

    let delta = nn_param_list
        .get("max_iter")
        .and_then(|v| v.as_real())
        .unwrap_or(0.001) as f32;

    let rho = nn_param_list
        .get("max_iter")
        .and_then(|v| v.as_real())
        .unwrap_or(1.0) as f32;

    NearestNeighbourParams {
        dist_metric,
        n_trees,
        search_budget,
        m,
        ef_construction,
        ef_search,
        max_iter,
        delta,
        rho,
    }
}

/// Wrapper function into the UMAP implementation in `manifolds-rs`
///
/// ### Params
///
/// ### Returns
#[allow(clippy::too_many_arguments)]
pub fn bxv_umap(
    data: MatRef<f32>,
    k: usize,
    min_dist: f32,
    spread: f32,
    n_dim: usize,
    ann_type: String,
    optim: String,
    init: String,
    seed: usize,
    verbose: bool,
) -> Mat<f32> {
    let optim_params =
        OptimParams::from_min_dist_spread(min_dist, spread, 1.0, 1.0, 500, 5, None, None, None);

    let res = umap(
        data,
        n_dim,
        k,
        optim,
        ann_type,
        init,
        &UmapParams::default(),
        None,
        Some(optim_params),
        false,
        seed,
        verbose,
    );

    let ncol = res.len();
    let nrow = res[0].len();

    Mat::from_fn(nrow, ncol, |i, j| res[j][i])
}
