//! Rust <> R interface for spatial analysis methods.

use bixverse_rs::graph::spatial_graph::{make_weights_non_redundant, GraphCsr};
use bixverse_rs::prelude::*;
use bixverse_rs::spatial::sp_analysis::nhood_enrichment::{
    compute_nhood_enrichment, NhoodEnrichmentParams, NhoodEnrichmentRes,
};
use bixverse_rs::spatial::sp_validate::validate_adjacency;
use extendr_api::*;
use std::time::Instant;

use crate::spatial::utils::adjacency_from_r_lists;

////////////////////
// extendr Module //
////////////////////

extendr_module! {
    mod r_sp_analysis;
    fn rs_sp_nhood_enrichment;
}

////////////////////////
// Nhood enrichment //
////////////////////////

/// Neighbourhood enrichment between spot labels
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Asks whether each pair of labels shows up as spatial neighbours more or
/// less often than chance. The graph is fixed and only the label assignment is
/// permuted, so the null keeps the spatial structure and tests the labelling
/// against it.
///
/// Each undirected edge contributes once to a diagonal entry and to both
/// off-diagonal entries, which keeps the count matrix symmetric throughout.
///
/// @param indices List of integer vectors. Local 0-based neighbour indices,
/// one element per spot. The output of [bixverse::rs_sp_build_graph()].
/// @param weights List of numeric vectors. Edge weights, aligned element for
/// element with `indices`.
/// @param labels Integer. Label code per spot, 0-based and local, so
/// `length(labels)` must equal `length(indices)`. Codes must be in
/// `0:(length(label_levels) - 1)`.
/// @param label_levels Character. Name per code, in encoding order.
/// @param nhood_params List. With the following elements:
/// \itemize{
///   \item n_perm - Integer. Number of permutations.
///   \item symmetrise - Boolean. Average `Z` with its transpose at the end.
/// }
/// @param seed Integer. Seed for the permutations.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @return A list with the following elements
/// \itemize{
///   \item label_levels - Character. Levels in encoding order, matching the
///   row and column order of every matrix below.
///   \item observed - Numeric matrix, K x K. Observed edge counts.
///   \item perm_mean - Numeric matrix, K x K. Permutation mean.
///   \item perm_sd - Numeric matrix, K x K. Permutation standard deviation.
///   \item z_scores - Numeric matrix, K x K. Z-scores.
/// }
///
/// @export
///
/// @keywords internal
#[extendr]
#[allow(clippy::too_many_arguments)]
fn rs_sp_nhood_enrichment(
    indices: List,
    weights: List,
    labels: Vec<i32>,
    label_levels: Vec<String>,
    nhood_params: List,
    seed: usize,
    verbose: usize,
) -> extendr_api::Result<List> {
    let verbosity = parse_verbosity_level(verbose);
    let params = NhoodEnrichmentParams::from_r_list(nhood_params)?;

    let (neighbours, weights) = adjacency_from_r_lists(&indices, &weights)?;

    if labels.len() != neighbours.len() {
        return Err(Error::Other(format!(
            "`labels` ({}) must be as long as `indices` ({}).",
            labels.len(),
            neighbours.len()
        )));
    }

    let labels: Vec<u32> = labels.iter().map(|&x| x as u32).collect();

    // The crate takes the symmetric CSR rather than the adjacency pair, so the
    // collapse-and-scatter that `build_spatial_graph_csr` does internally
    // happens here instead. That bypasses `MoransI::new`, which is where the
    // adjacency would otherwise be validated, so do it explicitly: a 1-based or
    // out-of-range neighbour list panics out of bounds inside the CSR scatter.
    validate_adjacency(&neighbours, &weights, labels.len()).to_extendr()?;

    let non_redundant = make_weights_non_redundant(&neighbours, &weights);
    let graph = GraphCsr::from_non_redundant(&neighbours, &non_redundant);

    let start = Instant::now();

    let res: NhoodEnrichmentRes =
        compute_nhood_enrichment(&graph, &labels, &label_levels, &params, seed).to_extendr()?;

    if verbosity.normal_verbosity() {
        println!(
            "Neighbourhood enrichment: {} labels, {} permutations in {:.2?}",
            res.n_labels,
            params.n_perm,
            start.elapsed()
        );
    }

    Ok(list!(
        label_levels = res.label_levels,
        observed = faer_to_r_matrix(res.observed.as_ref()),
        perm_mean = faer_to_r_matrix(res.perm_mean.as_ref()),
        perm_sd = faer_to_r_matrix(res.perm_sd.as_ref()),
        z_scores = faer_to_r_matrix(res.z_scores.as_ref())
    ))
}
