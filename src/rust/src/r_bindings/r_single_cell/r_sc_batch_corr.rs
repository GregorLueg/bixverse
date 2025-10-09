use extendr_api::*;

use crate::single_cell::batch_corrections::*;

/// Calculate kBET type scores
///
/// @export
#[extendr]
fn rs_kbet(knn_mat: RMatrix<i32>, batch_vector: Vec<i32>, threshold: f64) -> Vec<bool> {
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

    kbet(&knn_matrix, &batches, threshold)
}

extendr_module! {
    mod r_sc_batch_corr;
    fn rs_kbet;
}
