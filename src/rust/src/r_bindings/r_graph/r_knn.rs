use extendr_api::prelude::*;

use crate::core::graph::knn::KnnLabPropGraph;

/// kNN label propagation
///
/// @description
/// The function is a helper function to do kNN label propagation. This can
/// be useful for semi-supervised tasks. It implements the label spreading
/// method.
///
/// @param edge_list Integer vector. In form of node_1, node_2, node_3, ...
/// which indicates alternating pairs (node_1, node_2), etc in terms of edges
/// @param one_hot_encoding Integer matrix. Each row represents a sample, the
/// columns the one-hot encodings. Everything 0 denotes the unlabelled data.
/// @param label_mask Boolean vector. Which of the samples do not have a label.
/// Needs to be same length as `nrow(one_hot_encoding)`.
/// @param alpha Numeric. Parameter that controls the spreading. Usually between
/// 0.9 to 0.95. Larger values drive further labelling, smaller values are more
/// conversative.
/// @param iterations For how many (max) iterations to run the algorithm.
/// @param tolerance If the value below this is reached, an early stop is
/// initialised
///
/// @return The matrix with the probabilities of being of a certain class
///
/// @export
#[extendr]
fn rs_knn_label_propagation(
    edge_list: &[i32],
    one_hot_encoding: RMatrix<i32>,
    label_mask: &[Rbool],
    alpha: f64,
    iterations: usize,
    tolerance: f64,
) -> RMatrix<f64> {
    let edge_list = edge_list
        .iter()
        .map(|x| (*x - 1) as usize)
        .collect::<Vec<usize>>();
    let label_mask = label_mask
        .iter()
        .map(|r_obj| r_obj.to_bool())
        .collect::<Vec<bool>>();

    // transform the one hot encoding
    let nrow = one_hot_encoding.nrows();
    let ncol = one_hot_encoding.ncols();
    let data = one_hot_encoding.data();

    let n_nodes = nrow;

    let one_hot_encoding: Vec<Vec<f32>> = (0..nrow)
        .map(|i| (0..ncol).map(|j| data[i + j * nrow] as f32).collect())
        .collect();

    let knn_prog = KnnLabPropGraph::from_edge_list(&edge_list, n_nodes);

    let new_labels = knn_prog.label_spreading(
        &one_hot_encoding,
        &label_mask,
        alpha as f32,
        iterations,
        tolerance as f32,
    );

    RMatrix::new_matrix(nrow, ncol, |r, c| new_labels[r][c] as f64)
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

extendr_module! {
    mod r_knn;
    fn rs_knn_label_propagation;
    fn rs_knn_mat_to_edge_list;
}
