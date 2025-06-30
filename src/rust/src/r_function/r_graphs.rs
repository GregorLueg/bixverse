use extendr_api::prelude::*;
use rayon::prelude::*;

use crate::helpers::graph::*;
use crate::helpers::linalg::{col_means, col_sds};
use crate::utils_rust::nested_vector_to_faer_mat;

/// Helper function to calculate permutation-based page rank scores
///
/// @param node_names String vector. Name of the graph nodes.
/// @param from String vector. The names of the `from` edges from the edge list.
/// @param to String vector. The names of the `to` edges from the edge list.
/// @param diffusion_scores List. The personalised vectors for the page rank reset
/// values. Each element must sum to 1 and be of same length of `node_names`!
/// @param undirected Boolean. Is this an undirected graph.
///
/// @return A list containing:
///  \itemize{
///   \item means - The mean personalised page-rank scores based on the permutations.
///   \item sd - The standard deviation of the personalised page-rank scores based on
///   permutations.
/// }
///
/// @export
#[extendr]
fn rs_page_rank_permutations(
    node_names: Vec<String>,
    from: Vec<String>,
    to: Vec<String>,
    diffusion_scores: List,
    undirected: bool,
) -> extendr_api::Result<List> {
    let graph = graph_from_strings(&node_names, &from, &to, undirected);

    let mut personalise_vecs: Vec<Vec<f64>> = Vec::with_capacity(diffusion_scores.len());

    for i in 0..diffusion_scores.len() {
        let data = diffusion_scores.elt(i).unwrap();
        let vector_data: Vec<f64> = data.as_real_vector().unwrap();
        personalise_vecs.push(vector_data);
    }

    let page_rank_res: Vec<Vec<f64>> = personalise_vecs
        .par_iter()
        .map(|diff| personalized_page_rank(&graph, 0.85, diff, 1000, Some(1e-7)))
        .collect();

    let matrix_result = nested_vector_to_faer_mat(page_rank_res, false);

    let means = col_means(matrix_result.as_ref());
    let sds = col_sds(matrix_result.as_ref());

    Ok(list!(means = means, sd = sds))
}

/// Rust version of calcaluting the personalised page rank
///
/// @param node_names String vector. Name of the graph nodes.
/// @param from String vector. The names of the `from` edges from the edge list.
/// @param to String vector. The names of the `to` edges from the edge list.
/// @param personalised Numerical vector. The reset values. They must sum to 1 and
/// be of same length of `node_names`!
/// @param undirected Boolean. Is this an undirected graph.
///
/// @return The personalised page rank values.
///
/// @export
#[extendr]
fn rs_page_rank(
    node_names: Vec<String>,
    from: Vec<String>,
    to: Vec<String>,
    personalised: &[f64],
    undirected: bool,
) -> Vec<f64> {
    let graph = graph_from_strings(&node_names, &from, &to, undirected);

    personalized_page_rank(&graph, 0.85, personalised, 1000, Some(1e-7))
}

extendr_module! {
    mod r_graphs;
    fn rs_page_rank_permutations;
    fn rs_page_rank;
}
