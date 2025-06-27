use extendr_api::prelude::*;

use crate::helpers_sparse::SparseColumnMatrix;
use crate::utils_r_rust::sparse_matrix_to_list;

/// @export
#[extendr]
fn rs_upper_triangle_to_sparse(cor_vector: &[f64], shift: usize, n: usize) -> List {
    let include_diagonal = shift != 1;
    let sparse = SparseColumnMatrix::from_upper_triangle_sym(cor_vector, n, include_diagonal);

    sparse_matrix_to_list(sparse)
}

extendr_module! {
    mod fun_struct;
    fn rs_upper_triangle_to_sparse;
}
