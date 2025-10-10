use rand::prelude::*;
use rustc_hash::{FxBuildHasher, FxHashSet};
use std::cmp::PartialOrd;
use std::fmt::Debug;
use std::hash::Hash;
use std::ops::AddAssign;

use faer::{concat, Mat, MatRef};
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

/// Get the maximum value of an array
///
/// ### Params
///
/// * `arr` - The array of values
///
/// ### Returns
///
/// The maximum value found in the array
pub fn array_max<T: PartialOrd + Copy>(arr: &[T]) -> T {
    let mut max_val = arr[0];
    for number in arr {
        if *number > max_val {
            max_val = *number
        }
    }
    max_val
}

/// Get the minimum value of an array
///
/// ### Params
///
/// * `arr` - The array of values
///
/// ### Returns
///
/// The minimum value found in the array
pub fn array_min<T: PartialOrd + Copy>(arr: &[T]) -> T {
    let mut min_val = arr[0];
    for number in arr {
        if *number < min_val {
            min_val = *number
        }
    }
    min_val
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

/// Standard deviation
///
/// ### Params
///
/// * `x` Slice of `f64`
///
/// ### Returns
///
/// The standard deviation
pub fn standard_deviation(x: &[f64]) -> f64 {
    let n = x.len() as f64;
    let mean: f64 = x.iter().sum::<f64>() / n;
    let variance = x.iter().map(|&val| (val - mean).powi(2)).sum::<f64>() / (n - 1.0);
    variance.sqrt()
}

/// String slice to FxHashSet
///
/// ### Params
///
/// * `x` - The string slice.
///
/// ### Returns
///
/// A HashSet with borrowed String values
pub fn string_vec_to_set(x: &[String]) -> FxHashSet<&String> {
    let mut set = FxHashSet::with_capacity_and_hasher(x.len(), FxBuildHasher);
    for s in x {
        set.insert(s);
    }
    set
}

/// Get unique elements from a slice of any hashable, equatable numeric type.
///
/// ### Params
///
/// * `vec` - The slice of numerical values.
///
/// ### Returns
///
/// The unique elements of `vec` as a Vec.
pub fn unique<T>(vec: &[T]) -> Vec<T>
where
    T: Copy + Eq + Hash + Debug,
{
    let mut set = FxHashSet::default();
    vec.iter()
        .filter(|&&item| set.insert(item))
        .cloned()
        .collect()
}

/// Calculate the cumulative sum over a vector
///
/// ### Params
///
/// * `x` - The slice of numerical values
///
/// ### Returns
///
/// The cumulative sum over the vector.
pub fn cumsum<T>(x: &[T]) -> Vec<T>
where
    T: Copy + Default + AddAssign<T>,
{
    let mut sum = T::default();
    x.iter()
        .map(|&val| {
            sum += val;
            sum
        })
        .collect()
}

/// Split a vector randomly into two chunks
///
/// Splits a vector randomly into two of [..x] and the other [x..]
///
/// ### Params
///
/// * `vec` - Slice of the vector you want to split
/// * `x` - Length of the first vector; the rest will be put into the second vector
/// * `seed` - Seed for reproducibility
///
/// ### Returns
///
/// A tuple of the pieces of the vector
pub fn split_vector_randomly(vec: &[f64], x: usize, seed: u64) -> (Vec<f64>, Vec<f64>) {
    let mut rng = StdRng::seed_from_u64(seed);
    let mut shuffled = vec.to_vec();
    shuffled.shuffle(&mut rng);

    let (first_set, second_set) = shuffled.split_at(x);

    (first_set.to_vec(), second_set.to_vec())
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

/// Transform a nested vector into a faer matrix
///
/// ### Params
///
/// * `nested_vec` - The nested vector
/// * `col_wise` - If set to `True` it will column bind (outer vector represents)
///   the columns. If set to `False` it will row bind (outer vector represents
///   the rows).
///
/// ### Returns
///
/// The row or column bound matrix.
pub fn nested_vector_to_faer_mat(nested_vec: Vec<Vec<f64>>, col_wise: bool) -> Mat<f64> {
    let (nrow, ncol) = if col_wise {
        (nested_vec[0].len(), nested_vec.len())
    } else {
        (nested_vec.len(), nested_vec[0].len())
    };

    let data = flatten_vector(nested_vec);

    if col_wise {
        Mat::from_fn(nrow, ncol, |i, j| data[i + j * nrow])
    } else {
        Mat::from_fn(nrow, ncol, |i, j| data[j + i * ncol])
    }
}

/// Create a diagonal matrix from vector values
///
/// ### Params
///
/// * `vec` - The vector of values to put in the diagonal of the matrix
///
/// ### Returns
///
/// The diagonal matrix as a faer Matrix.
pub fn faer_diagonal_from_vec(vec: Vec<f64>) -> Mat<f64> {
    let len = vec.len();
    Mat::from_fn(len, len, |row, col| if row == col { vec[row] } else { 0.0 })
}

/// Get the index positions of the upper triangle of a symmetric matrix
///
/// Function will panic if offset > 1.
///
/// ### Params
///
/// * `n_dim` - The dimensions of the symmetric matrix
/// * `offset` - Do you want to include the diagonal values (offset = 0) or exclude
///   them (offset = 1).
///
/// ### Returns
///
/// A tuple of the row and column index positions of the upper triangle of the
/// matirx
pub fn upper_triangle_indices(n_dim: usize, offset: usize) -> (Vec<usize>, Vec<usize>) {
    if offset >= n_dim {
        return (Vec::new(), Vec::new());
    }
    assert!(offset <= 1, "The offset should be 0 or 1");

    // Precise calculation of total elements
    let total_elements: usize = (0..n_dim)
        .map(|row| n_dim.saturating_sub(row + offset))
        .sum();

    let mut row_indices = Vec::with_capacity(total_elements);
    let mut col_indices = Vec::with_capacity(total_elements);

    for row in 0..n_dim {
        let start_col = row + offset;
        if start_col < n_dim {
            let end_col = n_dim;
            let elements_in_row = end_col - start_col;
            // Use repeat_n for better performance and clarity
            row_indices.extend(std::iter::repeat_n(row, elements_in_row));
            col_indices.extend(start_col..end_col);
        }
    }

    (row_indices, col_indices)
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

/// Slice out a single row and return the remaining matrix
///
/// ### Params
///
/// * `x` - The matrix from which to remove a single row
/// * `idx_to_remove` - The index of the row to remove.
///
/// ### Returns
///
/// The matrix minus the specified row.
pub fn mat_rm_row(x: MatRef<f64>, idx_to_remove: usize) -> Mat<f64> {
    assert!(
        idx_to_remove <= x.nrows(),
        "The specified index is larger than the matrix"
    );

    let total_rows = x.nrows();

    let res = if idx_to_remove == 0 {
        x.subrows(1, total_rows - 1).to_owned()
    } else if idx_to_remove == total_rows - 1 {
        x.subrows(0, total_rows - 1).to_owned()
    } else {
        let upper = x.subrows(0, idx_to_remove);
        let lower = x.subrows(idx_to_remove + 1, total_rows - idx_to_remove - 1);
        concat![[upper], [lower]]
    };

    res
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

/// Colbind a vector of faer Matrices
///
/// The function will panic if the number of rows of the matrices differ in
/// the vector
///
/// ### Params
///
/// * `matrices` - Vector of faer matrix to column bind
///
/// ### Returns
///
/// One column bound matrix from the initial matrices
pub fn colbind_matrices(matrices: &[Mat<f64>]) -> Mat<f64> {
    let nrows = matrices[0].nrows();
    let total_col = matrices.iter().map(|m| m.ncols()).sum();
    let mut result: Mat<f64> = Mat::zeros(nrows, total_col);
    let mut col_offset = 0;
    for matrix in matrices {
        assert_eq!(
            matrix.nrows(),
            nrows,
            "All matrices must have the same number of columns"
        );
        let ncols = matrix.ncols();
        for i in 0..nrows {
            for j in 0..ncols {
                result[(i, col_offset + j)] = matrix[(i, j)]
            }
        }
        col_offset += ncols;
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use faer::mat;

    fn create_test_matrix_5x5() -> Mat<f64> {
        Mat::from_fn(5, 5, |i, j| (i * 5 + j + 1) as f64)
    }

    #[test]
    fn test_basic_slice_creation() {
        let matrix = create_test_matrix_5x5();
        let matrix_ref = matrix.as_ref();

        let row_indices = [0, 2, 4];
        let col_indices = [0, 2, 4];

        let slice_view = MatSliceView::new(matrix_ref, &row_indices, &col_indices);

        assert_eq!(slice_view.nrows(), 3);
        assert_eq!(slice_view.ncols(), 3);
    }

    #[test]
    fn test_slice_values() {
        let matrix = create_test_matrix_5x5();
        let matrix_ref = matrix.as_ref();

        // Original matrix:
        // [ 1.0,  2.0,  3.0,  4.0,  5.0]
        // [ 6.0,  7.0,  8.0,  9.0, 10.0]
        // [11.0, 12.0, 13.0, 14.0, 15.0]
        // [16.0, 17.0, 18.0, 19.0, 20.0]
        // [21.0, 22.0, 23.0, 24.0, 25.0]

        let row_indices = [0, 2, 4];
        let col_indices = [0, 2, 4];

        let slice_view = MatSliceView::new(matrix_ref, &row_indices, &col_indices);
        let sliced_matrix = slice_view.to_owned();

        // Expected result:
        // [ 1.0,  3.0,  5.0]
        // [11.0, 13.0, 15.0]
        // [21.0, 23.0, 25.0]

        let expected_result = mat![[1.0, 3.0, 5.0], [11.0, 13.0, 15.0], [21.0, 23.0, 25.0]];

        assert_eq!(sliced_matrix, expected_result);
    }

    #[test]
    fn test_inverted_order() {
        let matrix = create_test_matrix_5x5();
        let matrix_ref = matrix.as_ref();

        // Invert the order of indices
        let row_indices = [4, 2, 0];
        let col_indices = [4, 2, 0];

        let slice_view = MatSliceView::new(matrix_ref, &row_indices, &col_indices);
        let sliced_matrix = slice_view.to_owned();

        // Create expected matrix (inverted):
        // [25.0, 23.0, 21.0]
        // [15.0, 13.0, 11.0]
        // [ 5.0,  3.0,  1.0]
        let expected = mat![[25.0, 23.0, 21.0], [15.0, 13.0, 11.0], [5.0, 3.0, 1.0]];

        assert_eq!(sliced_matrix, expected);
    }

    #[test]
    fn test_duplicated_indices() {
        let matrix = create_test_matrix_5x5();
        let matrix_ref = matrix.as_ref();

        let row_indices: Vec<usize> = vec![0, 0, 2];
        let col_indices: Vec<usize> = vec![1, 1, 3];

        let slice_view = MatSliceView::new(matrix_ref, &row_indices, &col_indices);
        let sliced_matrix = slice_view.to_owned();

        // Create expected matrix (with duplicated indices):
        // [ 2.0,  2.0,  4.0]
        // [ 2.0,  2.0,  4.0]
        // [12.0, 12.0, 14.0]
        let expected = mat![[2.0, 2.0, 4.0], [2.0, 2.0, 4.0], [12.0, 12.0, 14.0]];

        assert_eq!(sliced_matrix, expected);
    }
}
