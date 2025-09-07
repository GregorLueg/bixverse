use faer::traits::ComplexField;
use faer::{Mat, MatRef};

/////////////
// Helpers //
/////////////

/// Counts the zeroes in a given faer matrix
///
/// ### Params
///
/// * `mat` - The respective faer matrix
///
/// ### Returns
///
/// A tuple with the first being the total zeroes, the second the zeroes per
/// row and the last element being the column zeroes.
pub fn count_zeroes<T>(mat: &MatRef<T>) -> (usize, Vec<usize>, Vec<usize>)
where
    T: Default + PartialEq,
{
    let (nrow, ncol) = mat.shape();
    let mut total_zeroes = 0_usize;
    let mut row_zeroes = vec![0_usize; nrow];
    let mut col_zeroes = vec![0_usize; ncol];

    let zero = T::default();

    #[allow(clippy::needless_range_loop)]
    for j in 0..ncol {
        for i in 0..nrow {
            let val = unsafe { mat.get_unchecked(i, j) };
            if *val == zero {
                total_zeroes += 1;
                row_zeroes[i] += 1;
                col_zeroes[j] += 1;
            }
        }
    }

    (total_zeroes, row_zeroes, col_zeroes)
}

////////////////
// Structures //
////////////////

/// Structure for SparseColumnMatrices
///
/// ### Fields
///
/// * `data` - Vector with the data.
/// * `row_indices` - The row indices of the data.
/// * `col_ptrs` - The column pointers of the data.
/// * `ncol` - Original number of columns.
/// * `nrow` - Original number of rows.
#[derive(Debug, Clone)]
pub struct SparseColumnMatrix<T> {
    pub data: Vec<T>,
    pub row_indices: Vec<usize>,
    pub col_ptrs: Vec<usize>,
    pub ncol: usize,
    pub nrow: usize,
}

#[allow(dead_code)]
impl<T> SparseColumnMatrix<T>
where
    T: Clone + Default + PartialEq + ComplexField + From<f64>,
{
    /// Generate a new sparse column matrix from values pre-computed data
    ///
    /// ### Params
    ///
    /// * `data` - Slice of the data.
    /// * `row_indices` - Slice of the row indices of the data.
    /// * `col_ptrs` - Slice of the column pointers of the data.
    /// * `ncol` - Original number of columns.
    /// * `nrow` - Original number of rows.
    pub fn new(
        data: &[T],
        row_indices: &[usize],
        col_ptrs: &[usize],
        ncol: usize,
        nrow: usize,
    ) -> Self {
        Self {
            data: data.to_vec(),
            row_indices: row_indices.to_vec(),
            col_ptrs: col_ptrs.to_vec(),
            ncol,
            nrow,
        }
    }

    /// Convert a faer dense matrix to sparse column format
    ///
    /// ### Params
    ///
    /// * `dense` - The original dense matrix.
    pub fn from_dense_matrix(dense: MatRef<T>) -> Self {
        let ncol = dense.ncols();
        let nrow = dense.nrows();

        let mut values = Vec::new();
        let mut row_indices = Vec::new();
        let mut col_ptrs = Vec::with_capacity(ncol + 1);

        col_ptrs.push(0_usize);

        for col in 0..ncol {
            for row in 0..nrow {
                let value = dense.get(row, col).clone();
                if value != T::default() {
                    values.push(value);
                    row_indices.push(row);
                }
            }
            col_ptrs.push(values.len());
        }

        Self {
            data: values,
            row_indices,
            col_ptrs,
            ncol,
            nrow,
        }
    }

    /// To a dense faer matrix
    ///
    /// ### Returns
    ///
    /// Returns a dense faer matrix.
    pub fn to_dense_matrix(&self) -> Mat<T> {
        let mut dense = Mat::zeros(self.nrow, self.ncol);

        for col in 0..self.ncol {
            let start = self.col_ptrs[col];
            let end = self.col_ptrs[col + 1];

            for idx in start..end {
                let row = self.row_indices[idx];
                let value = &self.data[idx];
                *dense.get_mut(row, col) = value.clone();
            }
        }

        dense
    }

    /// Return the number of non-zero values
    ///
    /// ### Returns
    ///
    /// Return the total number of NNZ values in the data
    pub fn nnz(&self) -> usize {
        self.data.len()
    }

    /// Create a sparse matrix from a row major upper triangle stored value
    ///
    /// This is a helper function to transform potentially sparse symmetric matrices
    /// stored as upper-triangles into a sparse matrix format in Rust.
    ///
    /// ### Params
    ///
    /// * `upper_triangle` - Represents the values of the upper triangle in
    ///   row major formant
    /// * `n` - Original nrows and ncols.
    /// * `include_diagonal` - Are the diagonal values included.
    pub fn from_upper_triangle_sym(upper_triangle: &[T], n: usize, include_diagonal: bool) -> Self {
        let mut values = Vec::new();
        let mut row_indices = Vec::new();
        let mut col_ptrs = Vec::with_capacity(n + 1);

        col_ptrs.push(0);

        for col in 0..n {
            for row in 0..n {
                let value = if row == col && !include_diagonal {
                    T::from(1.0)
                } else if row < col {
                    // Upper triangle (row < col)
                    let offset = if include_diagonal {
                        row * n - row * (row + 1) / 2 + col
                    } else {
                        row * (n - 1) - row * (row + 1) / 2 + col - 1
                    };
                    upper_triangle[offset].clone()
                } else if row == col {
                    // Diagonal (when included)
                    let offset = row * n - row * (row + 1) / 2 + col;
                    upper_triangle[offset].clone()
                } else {
                    // Lower triangle (col < row) - symmetric
                    let offset = if include_diagonal {
                        col * n - col * (col + 1) / 2 + row
                    } else {
                        col * (n - 1) - col * (col + 1) / 2 + row - 1
                    };
                    upper_triangle[offset].clone()
                };

                if value != T::default() {
                    values.push(value);
                    row_indices.push(row);
                }
            }
            col_ptrs.push(values.len());
        }

        Self {
            data: values,
            row_indices,
            col_ptrs,
            ncol: n,
            nrow: n,
        }
    }
}

//////////////////////////////
// Sparse format conversion //
//////////////////////////////

/// A type alias representing effect size results
///
/// ### Fields
///
/// * `0` - The data in the CSR format
/// * `1` - The row pointers
/// * `2` - The column indices
/// * `3` - An optional second data layer
pub type CsrData<T, U = T> = (Vec<T>, Vec<usize>, Vec<usize>, Option<Vec<U>>);

/// A type alias representing effect size results
///
/// ### Fields
///
/// * `0` - The data in the CSR format
/// * `1` - The column pointers
/// * `2` - The row indices
/// * `3` - An optional second data layer
pub type CscData<T, U = T> = (Vec<T>, Vec<usize>, Vec<usize>, Option<Vec<U>>);

/// Transform CSC stored data into CSR stored data
///
/// This version does a full memory copy of the data
///
/// ### Params
///
/// * `data` - The data stored in CSC format.
/// * `row_ind` - The row indices from the CSC format.
/// * `col_ptr` - The column pointers from the CSC format.
/// * `nrows` - The number of rows in the data.
/// * `data2` - An optional second data layer.
///
/// ### Returns
///
/// The data in CSR format, i.e., `CsrData`
#[allow(dead_code)]
pub fn csc_to_csr<T: Clone + Default, U: Clone + Default>(
    data: &[T],
    row_ind: &[usize],
    col_ptr: &[usize],
    nrows: usize,
    data2: Option<&[U]>,
) -> CsrData<T, U> {
    let nnz = data.len();
    let mut row_ptr = vec![0; nrows + 1];

    for &r in row_ind {
        row_ptr[r + 1] += 1;
    }

    for i in 0..nrows {
        row_ptr[i + 1] += row_ptr[i];
    }

    let mut csr_data = vec![T::default(); nnz];
    let mut csr_data2 = data2.map(|_| vec![U::default(); nnz]);
    let mut csr_col_ind = vec![0; nnz];
    let mut next = row_ptr[..nrows].to_vec();

    for col in 0..(col_ptr.len() - 1) {
        for idx in col_ptr[col]..col_ptr[col + 1] {
            let row = row_ind[idx];
            let pos = next[row];
            csr_data[pos] = data[idx].clone();
            csr_col_ind[pos] = col;

            // Add the second layer
            if let (Some(d2), Some(ref mut csr_d2)) = (data2, &mut csr_data2) {
                csr_d2[pos] = d2[idx].clone();
            }

            next[row] += 1;
        }
    }

    (csr_data, row_ptr, csr_col_ind, csr_data2)
}

/// Transform CSR stored data into CSC stored data
///
/// This version does a full memory copy of the data.
///
/// ### Params
///
/// * `data` - The data stored in CSR format.
/// * `col_ind` - The column indices from the CSR format.
/// * `row_ptr` - The row pointers from the CSR format.
/// * `ncols` - The number of columns in the data.
/// * `data2` - An optional second data layer.
///
/// ### Returns
///
/// The data in CSC format, i.e., `CscData`
pub fn csr_to_csc<T: Clone + Default, U: Clone + Default>(
    data: &[T],
    col_ind: &[usize],
    row_ptr: &[usize],
    ncols: usize,
    data2: Option<&[U]>,
) -> CscData<T, U> {
    let nnz = data.len();
    let mut col_ptr = vec![0; ncols + 1];

    // Count occurrences per column
    for &c in col_ind {
        col_ptr[c + 1] += 1;
    }

    // Cumulative sum to get column pointers
    for i in 0..ncols {
        col_ptr[i + 1] += col_ptr[i];
    }

    let mut csc_data = vec![T::default(); nnz];
    let mut csc_data2 = data2.map(|_| vec![U::default(); nnz]);
    let mut csc_row_ind = vec![0; nnz];
    let mut next = col_ptr[..ncols].to_vec();

    // Iterate through rows and place data in CSC format
    for row in 0..(row_ptr.len() - 1) {
        for idx in row_ptr[row]..row_ptr[row + 1] {
            let col = col_ind[idx];
            let pos = next[col];
            csc_data[pos] = data[idx].clone();
            csc_row_ind[pos] = row;

            // Add second layer if there
            if let (Some(d2), Some(ref mut csc_d2)) = (data2, &mut csc_data2) {
                csc_d2[pos] = d2[idx].clone();
            }

            next[col] += 1;
        }
    }

    (csc_data, col_ptr, csc_row_ind, csc_data2)
}

///////////
// Tests //
///////////

#[cfg(test)]
mod tests {
    use super::*;
    use faer::mat;

    #[test]
    fn test_dense_to_sparse_conversion() {
        let dense = mat![[1.0, 0.0, 3.0], [0.0, 2.0, 0.0], [4.0, 0.0, 5.0]];

        let sparse_obj = SparseColumnMatrix::from_dense_matrix(dense.as_ref());

        assert_eq!(sparse_obj.nrow, 3);
        assert_eq!(sparse_obj.ncol, 3);
        assert_eq!(sparse_obj.nnz(), 5);

        let expected_values = vec![1.0, 4.0, 2.0, 3.0, 5.0];

        assert_eq!(sparse_obj.data, expected_values);
    }

    #[test]
    fn test_dense_to_sparse_to_dense_conversion() {
        let dense = mat![[1.0, 0.0, 3.0], [0.0, 2.0, 0.0], [4.0, 0.0, 5.0]];

        let sparse_obj = SparseColumnMatrix::from_dense_matrix(dense.as_ref());

        let redense = sparse_obj.to_dense_matrix();

        assert_eq!(dense, redense);
    }

    #[test]
    fn test_raw_to_dense_conversion() {
        let data = vec![1.0, 4.0, 2.0, 3.0, 5.0];
        let row_indices: Vec<usize> = vec![0, 2, 1, 0, 2];
        let col_ptr: Vec<usize> = vec![0, 2, 3, 5];

        let sparse_obj = SparseColumnMatrix::new(&data, &row_indices, &col_ptr, 3, 3);

        let dense = mat![[1.0, 0.0, 3.0], [0.0, 2.0, 0.0], [4.0, 0.0, 5.0]];

        let redense = sparse_obj.to_dense_matrix();

        assert_eq!(dense, redense);
    }

    #[test]
    fn test_from_upper_triangle_vec_with_diagonal() {
        // Upper triangle with diagonal: [1, 0.8, 0.6, 1, 0.3, 1]
        // Represents correlation matrix:
        // [1.0, 0.8, 0.6]
        // [0.8, 1.0, 0.3]
        // [0.6, 0.3, 1.0]
        let upper_tri = vec![1.0, 0.8, 0.6, 1.0, 0.3, 1.0];
        let sparse = SparseColumnMatrix::from_upper_triangle_sym(&upper_tri, 3, true);

        let dense = sparse.to_dense_matrix();
        let expected = mat![[1.0, 0.8, 0.6], [0.8, 1.0, 0.3], [0.6, 0.3, 1.0]];

        assert_eq!(dense, expected);
    }

    #[test]
    fn test_from_upper_triangle_vec_without_diagonal() {
        // Upper triangle without diagonal: [0.8, 0.6, 0.3]
        // With implied diagonal of 1's, represents:
        // [1.0, 0.8, 0.6]
        // [0.8, 1.0, 0.3]
        // [0.6, 0.3, 1.0]
        let upper_tri = vec![0.8, 0.6, 0.3];
        let sparse = SparseColumnMatrix::from_upper_triangle_sym(&upper_tri, 3, false);

        let dense = sparse.to_dense_matrix();
        let expected = mat![[1.0, 0.8, 0.6], [0.8, 1.0, 0.3], [0.6, 0.3, 1.0]];

        assert_eq!(dense, expected);
    }

    #[test]
    fn test_csc_to_csr_conversion_single_layer() {
        // Test matrix:
        // [1 0 2]
        // [0 3 0]
        // [4 0 5]
        // CSC format
        let data = vec![1, 4, 3, 2, 5];
        let row_ind = vec![0, 2, 1, 0, 2];
        let col_ptr = vec![0, 2, 3, 5];
        let (csr_data, csr_row_ptr, csr_col_ind, csr_data2) =
            csc_to_csr(&data, &row_ind, &col_ptr, 3, None::<&[i32]>);

        // Expected CSR format
        assert_eq!(csr_data, vec![1, 2, 3, 4, 5]);
        assert_eq!(csr_col_ind, vec![0, 2, 1, 0, 2]);
        assert_eq!(csr_row_ptr, vec![0, 2, 3, 5]);
        assert_eq!(csr_data2, None);
    }

    #[test]
    fn test_csc_to_csr_conversion_dual_layer() {
        // Test matrix:
        // [1 0 2]    [10 0  20]
        // [0 3 0]    [0  30 0 ]
        // [4 0 5]    [40 0  50]
        // CSC format for both layers
        let data = vec![1, 4, 3, 2, 5];
        let data2 = vec![1.0, 4.0, 3.0, 2.0, 5.0];
        let row_ind = vec![0, 2, 1, 0, 2];
        let col_ptr = vec![0, 2, 3, 5];
        let (csr_data, csr_row_ptr, csr_col_ind, csr_data2) =
            csc_to_csr(&data, &row_ind, &col_ptr, 3, Some(&data2));

        // Expected CSR format
        assert_eq!(csr_data, vec![1, 2, 3, 4, 5]);
        assert_eq!(csr_col_ind, vec![0, 2, 1, 0, 2]);
        assert_eq!(csr_row_ptr, vec![0, 2, 3, 5]);
        assert_eq!(csr_data2, Some(vec![1.0, 2.0, 3.0, 4.0, 5.0]));
    }

    #[test]
    fn test_csr_to_csc_conversion_single_layer() {
        // Test matrix:
        // [1 0 2]
        // [0 3 0]
        // [4 0 5]
        // CSR format
        let data = vec![1, 2, 3, 4, 5];
        let col_ind = vec![0, 2, 1, 0, 2];
        let row_ptr = vec![0, 2, 3, 5];
        let (csc_data, csc_col_ptr, csc_row_ind, csc_data2) =
            csr_to_csc(&data, &col_ind, &row_ptr, 3, None::<&[i32]>);

        // Expected CSC format
        assert_eq!(csc_data, vec![1, 4, 3, 2, 5]);
        assert_eq!(csc_row_ind, vec![0, 2, 1, 0, 2]);
        assert_eq!(csc_col_ptr, vec![0, 2, 3, 5]);
        assert_eq!(csc_data2, None);
    }

    #[test]
    fn test_csr_to_csc_conversion_dual_layer() {
        // Test matrix:
        // [1 0 2]    [10 0  20]
        // [0 3 0]    [0  30 0 ]
        // [4 0 5]    [40 0  50]
        // CSR format for both layers
        let data = vec![1, 2, 3, 4, 5];
        let data2 = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let col_ind = vec![0, 2, 1, 0, 2];
        let row_ptr = vec![0, 2, 3, 5];
        let (csc_data, csc_col_ptr, csc_row_ind, csc_data2) =
            csr_to_csc(&data, &col_ind, &row_ptr, 3, Some(&data2));

        // Expected CSC format
        assert_eq!(csc_data, vec![1, 4, 3, 2, 5]);
        assert_eq!(csc_row_ind, vec![0, 2, 1, 0, 2]);
        assert_eq!(csc_col_ptr, vec![0, 2, 3, 5]);
        assert_eq!(csc_data2, Some(vec![1.0, 4.0, 3.0, 2.0, 5.0]));
    }
}
