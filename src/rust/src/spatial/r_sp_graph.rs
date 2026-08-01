//! Rust <> R interface for the spatial neighbourhood graph.

use bixverse_rs::graph::spatial_graph::{compute_node_degree, make_weights_non_redundant};
use bixverse_rs::prelude::*;
use bixverse_rs::spatial::sp_graph::{build_spatial_graph, SpatialAdjacency, SpatialGraphParams};
use extendr_api::*;
use std::time::Instant;

use crate::spatial::utils::{coords_from_r_matrix, neighbours_to_r_list, weights_to_r_list};

////////////////////
// extendr Module //
////////////////////

extendr_module! {
    mod r_sp_graph;
    fn rs_sp_build_graph;
}

///////////
// Graph //
///////////

/// Build the spatial neighbourhood graph for one sample
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Maps per-spot coordinates onto the `(neighbours, weights)` adjacency pair
/// that every other spatial statistic in this package consumes. Four layouts:
/// the two Visium lattices are exact table look-ups and never compute a
/// distance, kNN and fixed-radius cover irregular geometries like Xenium or
/// MERFISH.
///
/// The returned indices are **local**: entry `i` describes the spot at row `i`
/// of `coords`, and the indices inside it position against that same matrix.
/// They are not global spot ids.
///
/// @param coords Numeric matrix. Two columns, `x` then `y`, in full-resolution
/// pixels. One row per spot.
/// @param graph_params List. Layout and weighting. Flat list with:
/// \itemize{
///   \item layout - String. One of `c("hex", "square", "knn", "radius")`.
///   \item array_row, array_col - Integer. Lattice coordinates, needed by
///   `"hex"` and `"square"`.
///   \item connectivity - Integer. `4L` or `8L`, needed by `"square"`.
///   \item k - Integer. Neighbours per spot, needed by `"knn"`.
///   \item radius - Float. Cut-off in full-res pixels, needed by `"radius"`.
///   \item weighting - String. One of
///   `c("binary", "inverse_distance", "gaussian")`.
///   \item power - Float. Decay exponent for `"inverse_distance"`.
///   \item bandwidth - Float or `NULL`. Gaussian bandwidth; `NULL` resolves to
///   the median nearest-neighbour distance.
///   \item row_normalise - Boolean. spdep `style = "W"`.
/// }
/// @param knn_params List or `NULL`. ANN settings for the `"knn"` layout, see
/// [bixverse::params_sc_knn()]. Ignored by every other layout. `k` and the
/// distance metric are overridden from the layout and the 2D geometry.
/// @param seed Integer. Seed for the ANN index, for reproducibility.
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @return A list with the following elements
/// \itemize{
///   \item indices - List of integer vectors. Local 0-based neighbour indices,
///   sorted ascending, one element per spot.
///   \item weights - List of numeric vectors. Edge weights, aligned element
///   for element with `indices`.
///   \item degree - Numeric. Weighted degree per spot, taken over the
///   symmetrised graph `W + t(W)` that every statistic downstream consumes.
///   A reciprocal edge therefore contributes `w_ij + w_ji`, so a binary
///   layout with mutual neighbours gives twice the neighbour count.
/// }
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_sp_build_graph(
    coords: RMatrix<f64>,
    graph_params: List,
    knn_params: Nullable<List>,
    seed: usize,
    verbose: usize,
) -> extendr_api::Result<List> {
    let verbosity = parse_verbosity_level(verbose);
    let coordinates = coords_from_r_matrix(&coords)?;
    let params = SpatialGraphParams::from_r_list(graph_params)?;

    let knn_params = match knn_params {
        Nullable::NotNull(list) => Some(KnnParams::from_r_list(list)?),
        Nullable::Null => None,
    };

    let start = Instant::now();

    let (neighbours, weights): SpatialAdjacency =
        build_spatial_graph(&coordinates, &params, knn_params.as_ref(), seed).to_extendr()?;

    // `compute_node_degree` wants the non-redundant form, so the degree is the
    // row sum of the symmetric CSR the statistics actually consume. Handing it
    // the raw weights would double-count every reciprocal edge on one side.
    let degree = compute_node_degree(
        &neighbours,
        &make_weights_non_redundant(&neighbours, &weights),
    );

    if verbosity.normal_verbosity() {
        let edges: usize = neighbours.iter().map(|nb| nb.len()).sum();
        println!(
            "Spatial graph: {} spots, {} directed edges in {:.2?}",
            neighbours.len(),
            edges,
            start.elapsed()
        );
    }

    Ok(list!(
        indices = neighbours_to_r_list(&neighbours),
        weights = weights_to_r_list(&weights),
        degree = degree.iter().map(|&x| x as f64).collect::<Vec<f64>>()
    ))
}
