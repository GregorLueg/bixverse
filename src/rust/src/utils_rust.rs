use std::cmp::PartialOrd;
use std::collections::HashSet;
use std::fmt::Debug;
use std::hash::Hash;

use faer::Mat;

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

/// Transform a nested vector into a faer matrix
pub fn nested_vector_to_faer_mat(nested_vec: Vec<Vec<f64>>) -> faer::Mat<f64> {
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
