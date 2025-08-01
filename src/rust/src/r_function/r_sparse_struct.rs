use extendr_api::prelude::*;

use crate::helpers::structs_sparse::SparseColumnMatrix;
use crate::utils::r_rust_interface::sparse_matrix_to_list;

/// Generate sparse data from an upper triangle
///
/// @description This function takes the values from an upper triangle matrix
/// the shift and the nrows/ncols and returns a list.
///
/// @param value Numeric vector. The upper triangle values.
/// @param shift Integer Did you apply a shift to remove the diagonal values?
/// @param n Integer. The number of columns/rows in the symmetric matrix.
///
/// @return A list containing:
///  \itemize{
///   \item data - A vector of lists with the elements. (Related to the way
///   Robj are stored in Rust.)
///   \item row_indices - A vector of integers with the row indices.
///   \item col_ptr - A vector of integers with the column pointers.
/// }
///
/// @export
#[extendr]
fn rs_upper_triangle_to_sparse(value: &[f64], shift: usize, n: usize) -> List {
    let include_diagonal = shift != 1;
    let sparse = SparseColumnMatrix::from_upper_triangle_sym(value, n, include_diagonal);

    sparse_matrix_to_list(sparse)
}

extendr_module! {
    mod r_sparse_struct;
    fn rs_upper_triangle_to_sparse;
}
