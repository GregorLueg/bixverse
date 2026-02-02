use std::cmp::PartialOrd;
use std::fmt::Debug;

use faer::{Mat, MatRef};
use faer_entity::Entity;

//////////////////
// VECTOR STUFF //
//////////////////

/// Flatten a nested vector
///
/// ### Params
///
/// * `vec` - The vector to flatten
///
/// ### Returns
///
/// The flattened vector
pub fn flatten_vector<I, T>(vec: I) -> Vec<T>
where
    I: IntoIterator,
    I::Item: IntoIterator<Item = T>,
{
    vec.into_iter().flatten().collect()
}

/// Get the maximum and minimum value of an array
///
/// ### Params
///
/// * `arr` - The array of values
///
/// ### Returns
///
/// Tuple of values with the first being the minimum and the second the maximum
pub fn array_max_min<T: PartialOrd + Copy>(arr: &[T]) -> (T, T) {
    let mut min_val = arr[0];
    let mut max_val = arr[0];
    for number in arr {
        if *number < min_val {
            min_val = *number
        }
        if *number > max_val {
            max_val = *number
        }
    }

    (min_val, max_val)
}

//////////////////
// MATRIX STUFF //
//////////////////

/// Matrix slice view
///
/// Structure to help creating sub slices of a given matrix when needed.
/// Due to memory structure, slicing creates deep copies, but this function
/// avoids generating them until the last possible moment.
///
/// ### Fields
///
/// * `data` - The faer MatRef (original matrix)
/// * `row_indices` - The row indices you want to slice out.
/// * `col_indices` - The col indices you want to slice out.
#[derive(Clone, Debug)]
#[allow(dead_code)]
pub struct MatSliceView<'a, 'r, 'c, E: Entity> {
    data: MatRef<'a, E>,
    row_indices: &'r [usize],
    col_indices: &'c [usize],
}

#[allow(dead_code)]
impl<'a, 'r, 'c, E: Entity> MatSliceView<'a, 'r, 'c, E> {
    /// Generate a new MatSliceView
    ///
    /// This function will panic if you try to select indices larger than the
    /// underlying matrix.
    ///
    /// ### Params
    ///
    /// * `data` - The original MatRef from which you want to slice out data
    /// * `row_indices` - The row indices you want to slice out.
    /// * `col_indices` - The col indices you want to slice out.
    pub fn new(data: MatRef<'a, E>, row_indices: &'r [usize], col_indices: &'c [usize]) -> Self {
        let max_col_index = col_indices.iter().max().copied().unwrap_or(0);
        let max_row_index = row_indices.iter().max().copied().unwrap_or(0);

        assert!(
            max_col_index < data.ncols(),
            "You selected indices larger than ncol."
        );
        assert!(
            max_row_index < data.nrows(),
            "You selected indices larger than nrow."
        );

        Self {
            data,
            row_indices,
            col_indices,
        }
    }

    /// Return the number of rows
    ///
    /// ### Returns
    ///
    /// Number of rows
    pub fn nrows(&self) -> usize {
        self.row_indices.len()
    }

    /// Return the number of columns
    ///
    /// ### Returns
    ///
    /// Number of columns
    pub fn ncols(&self) -> usize {
        self.col_indices.len()
    }

    /// Return an owned matrix.
    ///
    /// Deep copying cannot be circumvented due to memory accessing at this point.
    /// Subsequent matrix algebra needs a continouos view into memory.
    ///
    /// ### Returns
    ///
    /// Owned sliced matrix for subsequent usage.
    pub fn to_owned(&self) -> Mat<E> {
        Mat::from_fn(self.nrows(), self.ncols(), |i, j| {
            *self.data.get(self.row_indices[i], self.col_indices[j])
        })
    }
}

/// Create from the upper triangle values a symmetric matrix
///
/// Generates the full dense matrix of values representing the upper triangle
/// of a symmetric matrix.
///
/// ### Params
///
/// * `data` - Slice of the values
/// * `shift` - Was the diagonal included (= 0) or not (= 1). If not included,
///   the diagonal is set to 1.
/// * `n` - Original dimension of the symmetric matrix.
///
/// ### Return
///
/// The symmetric, dense matrix.
pub fn upper_triangle_to_sym_faer(data: &[f64], shift: usize, n: usize) -> faer::Mat<f64> {
    let mut mat = Mat::<f64>::zeros(n, n);
    let mut idx = 0;
    for i in 0..n {
        for j in i..n {
            if shift == 1 && i == j {
                mat[(i, j)] = 1_f64
            } else {
                mat[(i, j)] = data[idx];
                mat[(j, i)] = data[idx];
                idx += 1;
            }
        }
    }

    mat
}

/// Store the upper triangle values as a flat vector from a faer matrix
///
/// ### Params
///
/// * `x` The faer matrix
/// * `shift` Shall the diagonal be included (shift = 0) or not (shift = 1).
///
/// ### Returns
///
/// A vector representing the upper triangle values (row major ordered)
pub fn faer_mat_to_upper_triangle(x: MatRef<f64>, shift: usize) -> Vec<f64> {
    assert!(shift <= 1, "The shift should be 0 or 1");

    let n = x.ncols();
    let total_elements = if shift == 0 {
        n * (n + 1) / 2
    } else {
        n * (n - 1) / 2
    };
    let mut vals = Vec::with_capacity(total_elements);
    for i in 0..n {
        let start_j = i + shift;
        for j in start_j..n {
            vals.push(*x.get(i, j));
        }
    }

    vals
}

/// Rowbind a vector of faer Matrices
///
/// The function will panic if the number of columns of the matrices differ in
/// the vector
///
/// ### Params
///
/// * `matrices` - Vector of faer matrix to row bind
///
/// ### Returns
///
/// One row bound matrix from the initial matrices
#[allow(dead_code)]
pub fn rowbind_matrices(matrices: Vec<Mat<f64>>) -> Mat<f64> {
    let ncols = matrices[0].ncols();
    let total_row = matrices.iter().map(|m| m.nrows()).sum();
    let mut result: Mat<f64> = Mat::zeros(total_row, ncols);
    let mut row_offset = 0;
    for matrix in matrices {
        assert_eq!(
            matrix.ncols(),
            ncols,
            "All matrices must have the same number of columns"
        );
        let nrows = matrix.nrows();
        for i in 0..nrows {
            for j in 0..ncols {
                result[(row_offset + i, j)] = matrix[(i, j)]
            }
        }
        row_offset += nrows;
    }

    result
}
