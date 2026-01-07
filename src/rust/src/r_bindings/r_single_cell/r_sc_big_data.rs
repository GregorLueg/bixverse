use extendr_api::*;
use faer::Mat;
use std::time::Instant;

use crate::single_cell::methods::big_data::*;
use crate::utils::r_rust_interface::{faer_to_r_matrix, r_matrix_to_faer_fp32};

extendr_module! {
    mod r_sc_big_data;
    fn rs_sc_knn_big_data;
}

/// @export
#[extendr]
fn rs_sc_knn_big_data(
    embd: RMatrix<f64>,
    knn_params: List,
    verbose: bool,
    seed: usize,
) -> extendr_api::RArray<i32, [usize; 2]> {
    let embd = r_matrix_to_faer_fp32(&embd);

    let knn_params = KnnParamsBigData::from_r_list(knn_params);

    let start_knn = Instant::now();

    let knn = generate_knn_quantised(embd.as_ref(), &knn_params, seed, verbose);

    let end_knn = start_knn.elapsed();

    if verbose {
        println!("KNN generation (for large data) done : {:.2?}", end_knn);
    }

    let index_mat = Mat::from_fn(embd.nrows(), knn_params.k, |i, j| knn[i][j] as i32);

    faer_to_r_matrix(index_mat.as_ref())
}
