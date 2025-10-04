use extendr_api::prelude::*;

use crate::core::base::cors_similarity::*;
use crate::core::base::rbf::rbf_gaussian_mat;
use crate::core::graph::graph_clustering::*;
use crate::utils::r_rust_interface::*;

/// Rust implementation of spectral clustering
///
/// @description Spectral clustering on a pre-calculated similarity matrix.
///
/// @param similarities Numerical matrix representing the similarities. Needs
/// to be symmetric!
/// @param k_neighbours Integer. Number of neighbours to consider in the kNN
/// graph generation
/// @param n_clusters Integer. Number of clusters to identify
/// @param max_iters Integer. Number of iterations for k-means clustering
/// @param seed Integer. Seed for reproducibility
///
/// @return A vector with the membership of the samples
///
/// @export
#[extendr]
fn rs_spectral_clustering_sim(
    similarities: RMatrix<f64>,
    k_neighbours: usize,
    n_clusters: usize,
    max_iters: usize,
    seed: usize,
) -> Vec<i32> {
    let sim = r_matrix_to_faer(&similarities);

    let membership = spectral_clustering(&sim, k_neighbours, n_clusters, max_iters, seed);

    membership
        .iter()
        .map(|x| (*x + 1_usize) as i32)
        .collect::<Vec<i32>>()
}

/// Rust implementation of spectral clustering
///
/// @param data Numerical matrix. The data to cluster. Rows = samples, columns =
/// features.
/// @param distance_type String. One of
/// `c("euclidean", "manhattan", "canberra", "cosine")`.
/// @param epsilon Numerical. The epsilon parameter for the Gaussian Radial
/// Basis function
/// @param k_neighbours Integer. Number of neighbours to consider in the kNN
/// graph generation
/// @param n_clusters Integer. Number of clusters to identify
/// @param max_iters Integer. Number of iterations for k-means clustering
/// @param seed Integer. Seed for reproducibility
///
/// @return A vector with the membership of the samples
///
/// @export
#[extendr]
fn rs_spectral_clustering(
    data: RMatrix<f64>,
    distance_type: String,
    epsilon: f64,
    k_neighbours: usize,
    n_clusters: usize,
    max_iters: usize,
    seed: usize,
) -> extendr_api::Result<Vec<i32>> {
    let data = &r_matrix_to_faer(&data).transpose();

    let dist_type = parse_distance_type(&distance_type)
        .ok_or_else(|| format!("Invalid Distance type: {}", distance_type))?;

    let res = match dist_type {
        DistanceType::L2Norm => column_pairwise_l2_norm(data),
        DistanceType::L1Norm => column_pairwise_l1_norm(data),
        DistanceType::Cosine => column_pairwise_cosine_dist(data),
        DistanceType::Canberra => column_pairwise_canberra_dist(data),
    };

    let sim = rbf_gaussian_mat(res.as_ref(), &epsilon);

    let membership = spectral_clustering(&sim.as_ref(), k_neighbours, n_clusters, max_iters, seed);

    let res = membership
        .iter()
        .map(|x| (*x + 1_usize) as i32)
        .collect::<Vec<i32>>();

    Ok(res)
}

extendr_module! {
    mod r_graph_clustering;
    fn rs_spectral_clustering_sim;
    fn rs_spectral_clustering;
}
