use bixverse_rs::graph::graph_label_propagations::*;
use extendr_api::prelude::*;
use std::collections::HashMap;

/////////////
// ExtendR //
/////////////

extendr_module! {
    mod r_knn;
    fn rs_knn_label_propagation;
    fn rs_knn_mat_to_edge_list;
    fn rs_knn_mat_to_edge_pairs;
}

/////////////
// Helpers //
/////////////

/// Helper for graph label propagation parameters
///
/// ### Params
///
/// * `alpha` - Parameter that controls the spreading
/// * `iter` - Number of iterations to run the algorithm for
/// * `tolerance` - Tolerance for early stopping
/// * `symmetrise` - Shall the graph be symmetrised
/// * `symmetry_strategy` - If weighted graph, which symmetrisation strategy
///   shall be used.
/// * `max_hops` - Optional usize. If provided, restricts label spreading to
///   nodes within this many hops of any labelled node. Nodes beyond this limit
///   are left as all-zeroes. If None, spreading is unrestricted.
struct LabelPropParams {
    alpha: f64,
    iter: usize,
    tolerance: f64,
    symmetrise: bool,
    symmetry_strategy: String,
    max_hops: Option<usize>,
}

impl LabelPropParams {
    fn from_list(r_list: List) -> Result<Self, extendr_api::Error> {
        let map: HashMap<&str, Robj> = r_list.try_into()?;

        let alpha = map.get("alpha").and_then(|v| v.as_real()).unwrap_or(0.9);

        let iter = map
            .get("iter")
            .and_then(|v| v.as_integer())
            .map(|v| v as usize)
            .unwrap_or(100);

        let tolerance = map
            .get("tolerance")
            .and_then(|v| v.as_real())
            .unwrap_or(1e-6);

        let symmetrise = map
            .get("symmetrise")
            .and_then(|v| v.as_logical())
            .map(|v| v.is_true())
            .unwrap_or(false);

        let symmetry_strategy = map
            .get("symmetry_strategy")
            .and_then(|v| v.as_str().map(|s| s.to_string()))
            .unwrap_or_else(|| "average".to_string());

        let max_hops = map
            .get("max_hops")
            .and_then(|v| v.as_integer())
            .map(|v| v as usize);

        Ok(Self {
            alpha,
            iter,
            tolerance,
            symmetrise,
            symmetry_strategy,
            max_hops,
        })
    }
}

////////////////////
// Main functions //
////////////////////

/// kNN label propagation
///
/// @description
/// The function is a helper function to do kNN label propagation. This can
/// be useful for semi-supervised tasks. It implements the label spreading
/// method.
///
/// @param from Integer vector. Source node indices for each edge.
/// @param to Integer vector. Target node indices for each edge. Must be the
/// same length as `from`.
/// @param one_hot_encoding Integer matrix. Each row represents a sample, the
/// columns the one-hot encodings. Everything 0 denotes the unlabelled data.
/// @param label_mask Boolean vector. Which of the samples do not have a label.
/// Needs to be same length as `nrow(one_hot_encoding)`.
/// @param weights Optional numeric vector. Edge weights for each pair in
/// `from`/`to`. Must have the same length as `from`. If NULL, all edges are
/// treated as unweighted.
/// @param label_prop_params List. Named list of parameters with the following
/// optional fields (defaults in parentheses):
/// \itemize{
///   \item \code{alpha} numeric, spreading strength (0.9)
///   \item \code{iter} integer, max iterations (100)
///   \item \code{tolerance} numeric, convergence threshold (1e-6)
///   \item \code{symmetrise} logical, symmetrise the graph (FALSE)
///   \item \code{symmetry_strategy} character, one of "average", "min", "max" ("average")
///   \item \code{max_hops} integer, restrict spreading radius (unrestricted)
/// }
///
/// @return The matrix with the probabilities of being of a certain class.
///
/// @export
#[extendr]
fn rs_knn_label_propagation(
    from: &[i32],
    to: &[i32],
    one_hot_encoding: RMatrix<i32>,
    label_mask: &[Rbool],
    weights: Nullable<Vec<f64>>,
    label_prop_params: List,
) -> Result<RMatrix<f64>, extendr_api::Error> {
    let params = LabelPropParams::from_list(label_prop_params)?;

    // deal with 1 indexing from R
    let from: Vec<usize> = from.iter().map(|x| (*x - 1) as usize).collect();
    let to: Vec<usize> = to.iter().map(|x| (*x - 1) as usize).collect();
    let label_mask: Vec<bool> = label_mask.iter().map(|r| r.to_bool()).collect();

    let nrow = one_hot_encoding.nrows();
    let ncol = one_hot_encoding.ncols();
    let data = one_hot_encoding.data();
    let n_nodes = nrow;

    let one_hot_encoding: Vec<Vec<f32>> = (0..nrow)
        .map(|i| (0..ncol).map(|j| data[i + j * nrow] as f32).collect())
        .collect();

    let knn_graph = match weights {
        Nullable::NotNull(w) => {
            let w: Vec<f32> = w.iter().map(|x| *x as f32).collect();
            let strategy = if params.symmetrise {
                parse_symmetry_strategy(&params.symmetry_strategy)
            } else {
                None
            };
            KnnLabPropGraph::from_weighted_node_pairs(&from, &to, &w, n_nodes, strategy)
        }
        Nullable::Null => KnnLabPropGraph::from_node_pairs(&from, &to, n_nodes, params.symmetrise),
    };

    let new_labels = knn_graph.label_spreading(
        &one_hot_encoding,
        &label_mask,
        params.alpha as f32,
        params.iter,
        params.tolerance as f32,
        params.max_hops,
    );

    Ok(RMatrix::new_matrix(nrow, ncol, |r, c| {
        new_labels[r][c] as f64
    }))
}

/// Flatten kNN matrix to edge list
///
/// @description
/// Helper function to leverage Rust to transform a kNN matrix into two vectors
/// of from, to
///
/// @param knn_mat Integer matrix. Rows represent the samples and the columns
/// the indices of the k-nearest neighbours.
/// @param one_index Boolean. If the original data is 0-index, shall 1-indexed
/// data be returned.
///
/// @return A flat vector representing the edge list.
///
/// @export
#[extendr]
fn rs_knn_mat_to_edge_list(knn_mat: RMatrix<i32>, one_index: bool) -> Vec<i32> {
    let samples = knn_mat.nrows();
    let neighbours = knn_mat.ncols();

    let mut res = Vec::with_capacity(samples * neighbours * 2);

    if one_index {
        for i in 0..samples {
            let from = (i + 1) as i32;
            for j in 0..neighbours {
                res.push(from);
                res.push(knn_mat[[i, j]] + 1);
            }
        }
    } else {
        for i in 0..samples {
            let from = i as i32;
            for j in 0..neighbours {
                res.push(from);
                res.push(knn_mat[[i, j]]);
            }
        }
    }

    res
}

/// Flatten kNN matrix to edge list
///
/// @description
/// Helper function to leverage Rust to transform a kNN matrix into an edge
/// list.
///
/// @param knn_mat Integer matrix. Rows represent the samples and the columns
/// the indices of the k-nearest neighbours.
/// @param one_index Boolean. If the original data is 0-index, shall 1-indexed
/// data be returned.
///
/// @return A list with the following elements
/// \itemize{
///   \item from - the from indices
///   \item to - the to indices
/// }
///
/// @export
#[extendr]
fn rs_knn_mat_to_edge_pairs(knn_mat: RMatrix<i32>, one_index: bool) -> List {
    let n = knn_mat.nrows();
    let k = knn_mat.ncols();

    let mut from: Vec<usize> = Vec::with_capacity(n * k);
    let mut to: Vec<usize> = Vec::with_capacity(n * k);

    for i in 0..n {
        for j in 0..k {
            let from_idx = if one_index { i + 1 } else { i };
            let to_idx = if one_index {
                (knn_mat[[i, j]] as usize) + 1
            } else {
                knn_mat[[i, j]] as usize
            };

            from.push(from_idx);
            to.push(to_idx);
        }
    }

    list![from = from, to = to]
}
