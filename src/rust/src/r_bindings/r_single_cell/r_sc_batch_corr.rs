use extendr_api::*;
use faer::Mat;

use crate::single_cell::bbknn::*;
use crate::utils::r_rust_interface::*;

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

/// BBKNN implementation in Rust
///
/// @description
/// This function implements the BBKNN algorithm from TO ADD
///
/// @param embd Numerical matrix. The embedding matrix to use to generate the
/// BBKNN parameters. Usually PCA. Rows represent cells.
/// @param batch_labels Integer vector. These represent to which batch a given
/// cell belongs.
/// @param bbknn_params List. Contains all of the BBKNN parameters.
/// @param seed Integer. Seed for reproducibility purposes.
/// @param verbose Boolean. Controls verbosity of the function.
///
/// @return A list of two lists representing the sparse matrix representation
/// of the distances and the connectivities.
///
/// @export
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

/// Reduce BBKNN matrix to Top X neighbours
///
/// @param indptr Integer vector. The index pointers of the underlying data.
/// @param indices Integer vector. The indices of the nearest neighbours.
/// @param cells_to_keep Integer. Number of nearest neighbours to keep.
///
/// @return A numerical matrix with the Top X neighbours per row. If
/// `cells_to_keep` is larger than the number of neighbours in the data, these
/// positions will be `NA`.
#[extendr]
fn rs_bbknn_filtering(
    indptr: Vec<i32>,
    indices: Vec<i32>,
    cells_to_keep: usize,
) -> extendr_api::RArray<f64, [usize; 2]> {
    let nrow = indptr.len() - 1;
    let ncol = cells_to_keep;
    let mut mat: Mat<f64> = Mat::from_fn(nrow, ncol, |_, _| f64::NAN);

    for i in 0..nrow {
        let start_i = indptr[i] as usize;
        let end_i = indptr[i + 1] as usize;
        let vals = end_i - start_i;
        let neighbours = if vals <= cells_to_keep {
            &indices[start_i..end_i]
        } else {
            &indices[start_i..(start_i + cells_to_keep)]
        };
        for (j, idx) in neighbours.iter().enumerate() {
            mat[(i, j)] = *idx as f64;
        }
    }
    faer_to_r_matrix(mat.as_ref())
}

extendr_module! {
    mod r_sc_batch_corr;
    fn rs_kbet;
    fn rs_bbknn;
    fn rs_bbknn_filtering;
}
