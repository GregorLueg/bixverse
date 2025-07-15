use extendr_api::prelude::*;

use rayon::prelude::*;
use std::sync::Arc;

use crate::helpers::graph::*;
use crate::helpers::linalg::{col_means, col_sds};
use crate::utils::general::nested_vector_to_faer_mat;

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

    // Arc version not needed here, as single run
    personalized_page_rank(&graph, 0.85, personalised, 1000, Some(1e-7))
}

/// Helper function to calculate permutation-based page rank scores
///
/// @description Helper function to calculate permutation based page-rank scores.
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
#[extendr]
fn rs_page_rank_permutations(
    node_names: Vec<String>,
    from: Vec<String>,
    to: Vec<String>,
    diffusion_scores: List,
    undirected: bool,
) -> extendr_api::Result<List> {
    let graph = graph_from_strings(&node_names, &from, &to, undirected);

    // Pre-process graph once
    let pagerank_graph = Arc::new(PageRankGraph::from_petgraph(&graph));

    let mut personalise_vecs: Vec<Vec<f64>> = Vec::with_capacity(diffusion_scores.len());
    for i in 0..diffusion_scores.len() {
        let data = diffusion_scores.elt(i)?;
        let vector_data: Vec<f64> = data.as_real_vector().unwrap();
        personalise_vecs.push(vector_data);
    }

    // Process in parallel with thread-local working memory
    let page_rank_res: Vec<Vec<f64>> = personalise_vecs
        .par_iter()
        .map_init(PageRankWorkingMemory::new, |working_memory, diff| {
            personalised_page_rank_optimised(
                &pagerank_graph,
                0.85,
                diff,
                1000,
                1e-7,
                working_memory,
            )
        })
        .collect();

    let matrix_result = nested_vector_to_faer_mat(page_rank_res, false);
    let means = col_means(matrix_result.as_ref());
    let sds = col_sds(matrix_result.as_ref());

    Ok(list!(means = means, sd = sds))
}

/// Helper function to calculate permutation-based tied page rank scores
///
/// @description Helper function to calculate permutation based tied page-rank
/// scores.
///
/// @param node_names String vector. Name of the graph nodes.
/// @param from String vector. The names of the `from` edges from the edge list.
/// @param to String vector. The names of the `to` edges from the edge list.
/// @param diffusion_scores_1 List. The first set of personalised vectors for
/// the page rank reset values. Each element must sum to 1 and be of same length
/// of `node_names`!
/// @param diffusion_scores_2 List. The second set of personalised vectors for
/// the page rank reset values. Each element must sum to 1 and be of same length
/// of `node_names`!
/// @param summarisation_fun String. One of `c("min", "max", "avg")`. Which type
/// of summarisation function to use to calculate the tied diffusion.
/// @param undirected Boolean. Is this an undirected graph.
///
/// @return A list containing:
///  \itemize{
///   \item means - The mean personalised page-rank scores based on the permutations.
///   \item sd - The standard deviation of the personalised page-rank scores based on
///   permutations.
/// }
#[extendr]
fn rs_page_rank_permutations_tied(
    node_names: Vec<String>,
    from: Vec<String>,
    to: Vec<String>,
    diffusion_scores_1: List,
    diffusion_scores_2: List,
    summarisation_fun: String,
    undirected: bool,
) -> extendr_api::Result<List> {
    assert!(
        diffusion_scores_1.len() == diffusion_scores_2.len(),
        "The two sets of random diffusion scores need to be the same length"
    );

    let graph_1 = graph_from_strings(&node_names, &from, &to, undirected);

    // For tied diffusion in the directed case, the directionality is reversed
    let graph_2 = if !undirected {
        graph_from_strings(&node_names, &to, &from, undirected)
    } else {
        graph_from_strings(&node_names, &from, &to, undirected)
    };

    // Pre-process graph once
    let pagerank_graph_1 = Arc::new(PageRankGraph::from_petgraph(&graph_1));
    let pagerank_graph_2 = Arc::new(PageRankGraph::from_petgraph(&graph_2));

    let mut personalise_vecs_1 = Vec::with_capacity(diffusion_scores_1.len());
    let mut personalise_vecs_2 = Vec::with_capacity(diffusion_scores_1.len());

    for i in 0..diffusion_scores_1.len() {
        personalise_vecs_1.push(diffusion_scores_1.elt(i)?.as_real_vector().unwrap());
        personalise_vecs_2.push(diffusion_scores_2.elt(i)?.as_real_vector().unwrap());
    }

    let summarisation_type: TiedSumType = parse_tied_sum(&summarisation_fun).unwrap();

    let tied_res: Vec<Vec<f64>> = personalise_vecs_1
        .into_par_iter()
        .zip(personalise_vecs_2.into_par_iter())
        .map_init(
            || (PageRankWorkingMemory::new(), PageRankWorkingMemory::new()),
            |(working_mem1, working_mem2), (diff1, diff2)| {
                let pr1 = personalised_page_rank_optimised(
                    &pagerank_graph_1,
                    0.85,
                    &diff1,
                    1000,
                    1e-7,
                    working_mem1,
                );
                let pr2 = personalised_page_rank_optimised(
                    &pagerank_graph_2,
                    0.85,
                    &diff2,
                    1000,
                    1e-7,
                    working_mem2,
                );

                pr1.iter()
                    .zip(pr2.iter())
                    .map(|(v1, v2)| match summarisation_type {
                        TiedSumType::Max => v1.max(*v2),
                        TiedSumType::Min => v1.min(*v2),
                        TiedSumType::Avg => (v1 + v2) * 0.5, // multiplication is faster
                    })
                    .collect()
            },
        )
        .collect();

    let matrix_result = nested_vector_to_faer_mat(tied_res, false);

    let means = col_means(matrix_result.as_ref());
    let sds = col_sds(matrix_result.as_ref());

    Ok(list!(means = means, sd = sds))
}

extendr_module! {
    mod r_graphs;
    fn rs_page_rank_permutations;
    fn rs_page_rank;
    fn rs_page_rank_permutations_tied;
}
