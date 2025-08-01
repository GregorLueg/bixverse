use faer::traits::ComplexField;
use faer::{Mat, MatRef};

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
    ///                      row major formant
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
}
