use extendr_api::*;

use crate::single_cell::batch_corrections::*;
use crate::utils::r_rust_interface::{r_matrix_to_faer_fp32, sparse_data_to_list};

/// Calculate kBET type scores
///
/// @description
/// The function takes in a kNN matrix and a batch vector indicating which
/// cell belongs to which batch. The function will check for the neighbourhood
/// of each cell if the proportion of represented batches are different from
/// the overall batch proportions. Good mixing of batches would mean very
/// cells have significant differences; bad mixing a lot of the batches
/// have bad mixing.
///
/// @param knn_mat Integer matrix. The rows represent the cells and the
/// columns the neighbour indices.
/// @param batch_vector Integer vector. The integers indicate to which
/// batch a given cell belongs.
///
/// @return A vector of p-values based on the ChiSquare statistic per cell.
///
/// @export
#[extendr]
fn rs_kbet(knn_mat: RMatrix<i32>, batch_vector: Vec<i32>) -> Vec<f64> {
    let n_cells = knn_mat.nrows();
    let k_neighbours = knn_mat.ncols();

    let knn_matrix: Vec<Vec<usize>> = (0..n_cells)
        .map(|i| {
            (0..k_neighbours)
                .map(|j| knn_mat[[i, j]] as usize)
                .collect()
        })
        .collect();

    // Convert batch_vector to Vec<usize>
    let batches: Vec<usize> = batch_vector.iter().map(|&x| x as usize).collect();

    kbet(&knn_matrix, &batches)
}

#[extendr]
fn rs_bbknn(
    embd: RMatrix<f64>,
    batch_labels: Vec<i32>,
    bbknn_params: List,
    seed: usize,
    verbose: bool,
) -> List {
    let bbknn_params = BbknnParams::from_r_list(bbknn_params);
    let embd = r_matrix_to_faer_fp32(&embd);
    let batch_labels = batch_labels
        .iter()
        .map(|x| *x as usize)
        .collect::<Vec<usize>>();

    let (distances, connectivities) =
        bbknn(embd.as_ref(), &batch_labels, &bbknn_params, seed, verbose);

    list!(
        distances = sparse_data_to_list(distances),
        connectivities = sparse_data_to_list(connectivities)
    )
}

extendr_module! {
    mod r_sc_batch_corr;
    fn rs_kbet;
    fn rs_bbknn;
}
