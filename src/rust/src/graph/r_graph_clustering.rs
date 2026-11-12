use bixverse_r::core::base::cors_similarity::*;
use bixverse_r::core::math::rbf::rbf_gaussian_mat;
use bixverse_r::graph::graph_clustering::*;
use bixverse_r::graph::graph_structures::edge_pairs_to_sparse_graph;
use bixverse_r::graph::leiden::{leiden_sparse_graph, LeidenParams};
use bixverse_r::ml::clustering::k_means::KMeansParamsWrappers;
use bixverse_r::prelude::*;
use extendr_api::prelude::*;

/////////////
// extendR //
/////////////

extendr_module! {
    mod r_graph_clustering;
    fn rs_spectral_clustering_sim;
    fn rs_spectral_clustering;
    fn rs_leiden_clustering;
}

///////////////
// Functions //
///////////////

/// Rust implementation of Leiden and Louvain clustering
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Modularity-based Leiden (or Louvain, with `refine = FALSE`) on an
/// undirected weighted graph. Mirrors
/// `igraph::cluster_leiden(objective_function = "modularity")`. The result
/// depends on the seed only, not on the number of threads.
///
/// @param from Integer vector. 1-based source node of each edge.
/// @param to Integer vector. 1-based target node of each edge. Same length as
/// `from`.
/// @param weights Numeric vector. Edge weights, same length as `from`.
/// @param n_nodes Integer. Number of nodes in the graph.
/// @param leiden_params List. Named list with the optional fields (defaults
/// in parentheses):
/// \itemize{
///   \item `resolution` numeric, modularity resolution (1.0)
///   \item `beta` numeric, randomness of the refinement merges (0.01)
///   \item `n_iterations` integer, full Leiden iterations (2)
///   \item `refine` logical, `TRUE` for Leiden, `FALSE` for Louvain (TRUE)
///   \item `seed` integer, random seed (42)
/// }
///
/// @returns Integer vector with the 1-based community of each node; `1` is the
/// largest community.
///
/// @export
#[extendr]
fn rs_leiden_clustering(
    from: &[i32],
    to: &[i32],
    weights: &[f64],
    n_nodes: usize,
    leiden_params: List,
) -> extendr_api::Result<Vec<i32>> {
    if from.len() != to.len() || from.len() != weights.len() {
        return Err(Error::Other(format!(
            "from ({}), to ({}) and weights ({}) must have the same length.",
            from.len(),
            to.len(),
            weights.len()
        )));
    }
    let params = LeidenParams::from_r_list(leiden_params)?;

    let mut edges = Vec::with_capacity(2 * from.len());
    for (&f, &t) in from.iter().zip(to) {
        for x in [f, t] {
            if x < 1 || x as usize > n_nodes {
                return Err(Error::Other(format!(
                    "Node index {} is outside 1..={}.",
                    x, n_nodes
                )));
            }
            edges.push((x - 1) as usize);
        }
    }

    let graph = edge_pairs_to_sparse_graph(&edges, weights, n_nodes);
    let membership = leiden_sparse_graph(&graph, Some(params)).to_extendr()?;

    Ok(membership.iter().map(|&x| (x + 1) as i32).collect())
}

/// Rust implementation of spectral clustering
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Spectral clustering on a pre-calculated similarity matrix.
///
/// @param similarities Numerical matrix representing the similarities. Needs
/// to be symmetric!
/// @param k_neighbours Integer. Number of neighbours to consider in the kNN
/// graph generation
/// @param n_clusters Integer. Number of clusters to identify
/// @param max_iters Integer. Number of iterations for k-means clustering
/// @param seed Integer. Seed for reproducibility
///
/// @returns A vector with the 1-based cluster membership of the samples.
///
/// @export
#[extendr]
fn rs_spectral_clustering_sim(
    similarities: RMatrix<f64>,
    k_neighbours: usize,
    n_clusters: usize,
    max_iters: usize,
    seed: usize,
) -> Result<Vec<i32>, extendr_api::Error> {
    let sim = r_matrix_to_faer(&similarities);

    let params = KMeansParamsWrappers::new(max_iters, None, None);

    let membership =
        spectral_clustering(&sim, k_neighbours, n_clusters, Some(params), seed).to_extendr()?;

    let res = membership
        .iter()
        .map(|x| (*x + 1_usize) as i32)
        .collect::<Vec<i32>>();

    Ok(res)
}

/// Rust implementation of spectral clustering
///
/// @description
/// `r lifecycle::badge("experimental")`
/// This version can take in data as is and will calculate the distance matrix
/// internally and convert it into similarities with a Gaussian RBF.
///
/// @param data Numerical matrix. The data to cluster. Rows = samples, columns =
/// features.
/// @param distance_type String. One of
/// `c("euclidean", "manhattan", "canberra", "cosine", "correlation")`. Unknown
/// strings raise an error.
/// @param epsilon Numerical. The epsilon parameter for the Gaussian Radial
/// Basis function
/// @param k_neighbours Integer. Number of neighbours to consider in the kNN
/// graph generation
/// @param n_clusters Integer. Number of clusters to identify
/// @param max_iters Integer. Number of iterations for k-means clustering
/// @param seed Integer. Seed for reproducibility
///
/// @returns A vector with the 1-based cluster membership of the samples.
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
        DistanceType::Correlation => {
            let cor = bixverse_r::core::base::cors_similarity::column_pairwise_cor(data, false);
            faer::Mat::from_fn(cor.nrows(), cor.ncols(), |i, j| 1.0 - cor[(i, j)])
        }
    };

    let params = KMeansParamsWrappers::new(max_iters, None, None);

    let sim = rbf_gaussian_mat(res.as_ref(), &epsilon);

    let membership =
        spectral_clustering(&sim.as_ref(), k_neighbours, n_clusters, Some(params), seed)
            .to_extendr()?;

    let res = membership
        .iter()
        .map(|x| (*x + 1_usize) as i32)
        .collect::<Vec<i32>>();

    Ok(res)
}
