use std::cmp::PartialOrd;
use std::collections::HashSet;
use std::fmt::Debug;
use std::hash::Hash;

use faer::{concat, Mat, MatRef};

//////////////////
// VECTOR STUFF //
//////////////////

/// Flatten a nested vector
pub fn flatten_vector<I, T>(vec: I) -> Vec<T>
where
    I: IntoIterator,
    I::Item: IntoIterator<Item = T>,
{
    vec.into_iter().flatten().collect()
}

/// Get the maximum value of an array
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
pub fn array_min<T: PartialOrd + Copy>(arr: &[T]) -> T {
    let mut min_val = arr[0];
    for number in arr {
        if *number < min_val {
            min_val = *number
        }
    }
    min_val
}

/// Get the maximum and minimum value. First element is minimum;
/// second one is maximum.
pub fn array_f64_max_min<T: PartialOrd + Copy>(arr: &[T]) -> (T, T) {
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

/// String vector to HashSet
pub fn string_vec_to_set(x: &[String]) -> HashSet<&String> {
    let mut set = HashSet::with_capacity(x.len());
    for s in x {
        set.insert(s);
    }
    set
}

/// Generate the rank of a vector with tie correction.
pub fn rank_vector(vec: &[f64]) -> Vec<f64> {
    let mut vec_index: Vec<(f64, usize)> = vec
        .iter()
        .copied()
        .enumerate()
        .map(|(i, x)| (x, i))
        .collect();

    vec_index.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    let mut ranks = vec![0.0; vec.len()];
    let mut i = 0;

    while i < vec_index.len() {
        let value = vec_index[i].0;
        let mut j = i + 1;

        // Tie correction
        while j < vec_index.len() && vec_index[j].0 == value {
            j += 1;
        }

        let rank = (i + j - 1) as f64 / 2.0 + 1.0;

        vec_index[i..j].iter().for_each(|&(_, original_index)| {
            ranks[original_index] = rank;
        });

        i = j;
    }

    ranks
}

/// Get unique elements from a slice of any hashable, equatable numeric type.
pub fn unique<T>(vec: &[T]) -> Vec<T>
where
    T: Copy + Eq + Hash + Debug,
{
    let mut set = HashSet::new();
    vec.iter()
        .filter(|&&item| set.insert(item))
        .cloned()
        .collect()
}

/// Calculate the cumulative sum over a vector
pub fn cumsum(values: &[f64]) -> Vec<f64> {
    let mut sum = 0.0;
    values
        .iter()
        .map(|&x| {
            sum += x;
            sum
        })
        .collect()
}

//////////////////
// MATRIX STUFF //
//////////////////

/// Structure to generate a slice of the data.
/// Key issue is that the data of the matrix is scattered in memory, thus,
/// scattered data access and deep copying WILL be needed at some point...
/// This structure will materialise the matrix ONLY when needed. Also, heavy
/// use of life times, so individual vectors/data can outlive the rest.
#[derive(Clone, Debug)]
pub struct MatSliceView<'a, 'r, 'c> {
    data: MatRef<'a, f64>,
    row_indices: &'r [usize],
    col_indices: &'c [usize],
}

impl<'a, 'r, 'c> MatSliceView<'a, 'r, 'c> {
    /// Generate a new MatSliceView
    pub fn new(data: MatRef<'a, f64>, row_indices: &'r [usize], col_indices: &'c [usize]) -> Self {
        Self {
            data,
            row_indices,
            col_indices,
        }
    }

    /// Return the number of rows
    pub fn nrows(&self) -> usize {
        self.row_indices.len()
    }

    /// Return the number of columns
    pub fn ncols(&self) -> usize {
        self.col_indices.len()
    }

    /// Return an owned matrix. Deep copying CANNOT be circumvented due to
    /// memory accessing at this point. Subsequent matrix algebra needs a continouos
    /// view into memory.
    pub fn to_owned(&self) -> Mat<f64> {
        Mat::from_fn(self.nrows(), self.ncols(), |i, j| {
            self.data[(self.row_indices[i], self.col_indices[j])]
        })
    }
}

/// Transform a nested vector into a faer matrix
pub fn nested_vector_to_faer_mat(nested_vec: Vec<Vec<f64>>) -> Mat<f64> {
    let nrow = nested_vec[0].len();
    let ncol = nested_vec.len();
    let data = flatten_vector(nested_vec);
    Mat::from_fn(nrow, ncol, |i, j| data[i + j * nrow])
}

/// Create a diagonal matrix with the vector values in the diagonal and the rest being 0's
pub fn faer_diagonal_from_vec(vec: Vec<f64>) -> Mat<f64> {
    let len = vec.len();
    Mat::from_fn(len, len, |row, col| if row == col { vec[row] } else { 0.0 })
}

/// Get the index positions of the upper triangle of a symmetric matrix
pub fn upper_triangle_indices(n_dim: usize, offset: usize) -> (Vec<usize>, Vec<usize>) {
    if offset >= n_dim {
        return (Vec::new(), Vec::new());
    }

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

            // Use extend with iterator for better performance
            row_indices.extend(std::iter::repeat(row).take(elements_in_row));
            col_indices.extend(start_col..end_col);
        }
    }

    (row_indices, col_indices)
}

/// Create from the upper triangle values for a symmetric matrix the full
/// dense faer matrix.
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

/// Slice out a single row and return the remaining
pub fn mat_row_rm_row(x: MatRef<f64>, idx_to_remove: usize) -> Mat<f64> {
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

// /// Rowbind a vector of faer Matrices, assuming same column length for all of
// /// them
// pub fn rowbind_matrices(
//   matrices: Vec<Mat<f64>>
// ) -> Mat<f64> {
//   let ncols = matrices[0].ncols();
//   let total_row = matrices
//     .iter()
//     .map(|m| m.nrows())
//     .sum();
//   let mut result: Mat<f64> = Mat::zeros(total_row, ncols);
//   let mut row_offset = 0;
//   for matrix in matrices{
//     assert_eq!(matrix.ncols(), ncols, "All matrices must have the same number of columns");
//     let nrows = matrix.nrows();
//     for i in 0..nrows {
//       for j in 0..ncols {
//         result[(row_offset + i, j)] = matrix[(i, j)]
//       }
//     }
//     row_offset += nrows;
//   }

//   result
// }

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
