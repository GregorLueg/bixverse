use bixverse_rs::core::math::sparse::count_zeroes;
use bixverse_rs::prelude::*;
use extendr_api::prelude::*;

/// Generate sparse data from an upper triangle
///
/// @description This function takes the values from an upper triangle matrix
/// the shift and the nrows/ncols and returns a list.
///
/// @param value Numeric vector. The upper triangle values.
/// @param shift Integer Did you apply a shift to remove the diagonal values?
/// @param n Integer. The number of columns/rows in the symmetric matrix.
/// @param cs_type String. One of `c("csr", "csc")`. Which type of list to
/// return.
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
fn rs_upper_triangle_to_sparse(
    value: &[f64],
    shift: usize,
    n: usize,
    cs_type: &str,
) -> extendr_api::Result<List> {
    let include_diagonal = shift != 1;
    let cs_type = parse_compressed_sparse_format(cs_type)
        .ok_or_else(|| extendr_api::Error::Other("Invalid cs_type".into()))?;
    let sparse = CompressedSparseData::from_upper_triangle_sym(value, n, include_diagonal, cs_type);
    Ok(sparse_data_to_list(sparse))
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
