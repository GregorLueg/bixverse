use extendr_api::prelude::*;

use crate::core::data::sparse_structures::*;
use crate::utils::r_rust_interface::{r_matrix_to_faer, sparse_matrix_to_list};

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

/// Helper to get zero stats from a given matrix
///
/// @description
/// Calculates in a single matrix pass the total number of zeroes, the row
/// zeroes and column zeroes.
///
/// @param x Numeric matrix. The matrix for which to calculate the total zeroes,
/// column and row zeroes.
///
/// @returns A list with:
/// \itemize{
///     \item total_zeroes - Total zeroes in the matrix.
///     \item row_zeroes - Vector with number of zeroes per row.
///     \item col_zeroes - Vector with number of zeroes per column.
/// }
///
/// @export
#[extendr]
fn rs_count_zeroes(x: RMatrix<f64>) -> List {
    let mat = r_matrix_to_faer(&x);
    let (total_zeroes, row_zeroes, col_zeroes) = count_zeroes(&mat);

    list!(
        total_zeroes = total_zeroes,
        row_zeroes = row_zeroes,
        col_zeroes = col_zeroes
    )
}

extendr_module! {
    mod r_sparse;
    fn rs_upper_triangle_to_sparse;
    fn rs_count_zeroes;
}
