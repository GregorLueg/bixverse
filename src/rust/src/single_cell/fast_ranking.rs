use half::f16;
use rayon::prelude::*;

use crate::utils::traits::F16;

/// Helper function to rank specifically `F16` type slices
///
/// ### Params
///
/// * `vec` - Slice of `F16`
///
/// ### Returns
///
/// The ranked values as an f16 vector.
#[allow(dead_code)]
fn rank_f16(vec: &[F16]) -> Vec<f16> {
    let n = vec.len();
    if n == 0 {
        return Vec::new();
    }

    let mut indexed_values: Vec<(F16, usize)> = vec
        .iter()
        .copied()
        .enumerate()
        .map(|(i, v)| (v, i))
        .collect();

    indexed_values
        .sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));

    let mut ranks: Vec<f16> = vec![f16::from_f64(0.0); n];
    let mut i = 0;
    while i < n {
        let current_value = indexed_values[i].0;
        let start = i;
        while i < n && indexed_values[i].0 == current_value {
            i += 1;
        }
        let avg_rank = f16::from_f32((start + i + 1) as f32) / f16::from_f32(2.0);
        for j in start..i {
            ranks[indexed_values[j].1] = avg_rank;
        }
    }
    ranks
}

/// Fast ranking of CSR-type data for single cell
///
/// The function takes in CSR-style data (rows = cells, columns = genes) and
/// generates ranked versions of the data.
///
/// ### Params
///
/// * `row_ptr` - The row pointer in the given CSR data
/// * `col_indices` - The col indices of the data
/// * `data` - The normalised count data which to rank
/// * `nrow` - Number of rows (cells)
/// * `ncol` - Number of columns (genes)
/// * `transpose` - This boolean indicates if the outer vector represents rows
///   (when set to `false`) or columns (when set to `true`). For example in the
///   case of Mann Whitney stats-based DGEs you would want this to be set to
///   `true` while for AUCell type statistics to `false`
///
/// ### Return
///
/// A `Vec<Vec<f16>>` that pending the transpose parameter represents the ranked
/// data in one or the other manner.
#[allow(dead_code)]
pub fn fast_csr_ranking(
    row_ptr: &[usize],
    col_indices: &[u16],
    data: &[F16],
    nrow: usize,
    ncol: usize,
    transpose: bool,
) -> Vec<Vec<f16>> {
    let row_ranked: Vec<Vec<f16>> = (0..nrow)
        .into_par_iter()
        .map(|row_idx| {
            let start = row_ptr[row_idx];
            let end = row_ptr[row_idx + 1];
            let num_nonzeros = end - start;
            let num_zeros = ncol - num_nonzeros;

            if num_nonzeros == 0 {
                let zero_rank = f16::from_f32((1.0 + ncol as f32) / 2.0);
                return vec![zero_rank; ncol];
            }

            if num_zeros == 0 {
                let row_data = &data[start..end];
                return rank_f16(row_data);
            }

            let row_data = &data[start..end];
            let row_cols = &col_indices[start..end];
            let nonzero_ranks = rank_f16(row_data);
            let zero_rank = f16::from_f32((1.0 + num_zeros as f32) / 2.0);
            let mut result = vec![zero_rank; ncol];

            for (i, &col) in row_cols.iter().enumerate() {
                result[col as usize] = nonzero_ranks[i] + f16::from_f32(num_zeros as f32);
            }

            result
        })
        .collect();

    if transpose {
        (0..ncol)
            .map(|col_idx| {
                (0..nrow)
                    .map(|row_idx| row_ranked[row_idx][col_idx])
                    .collect()
            })
            .collect()
    } else {
        row_ranked
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::traits::F16;
    use half::f16;

    // Helper to create F16 from f32
    fn f16_vec(values: &[f32]) -> Vec<F16> {
        values.iter().map(|&v| F16::from_f32(v)).collect()
    }

    // Helper to convert f16 vec to f32 for easier comparison
    fn to_f32_vec(values: &[f16]) -> Vec<f32> {
        values.iter().map(|&v| v.to_f32()).collect()
    }

    #[test]
    fn test_simple_row_ranking() {
        // Matrix
        // [1.0, 0.0, 3.0]
        // [0.0, 2.0, 0.0]
        let row_ptr = vec![0, 2, 3];
        let col_indices = vec![0, 2, 1];
        let data = f16_vec(&[1.0, 3.0, 2.0]);

        let result = fast_csr_ranking(&row_ptr, &col_indices, &data, 2, 3, false);

        // Row 0: [1.0, 0.0, 3.0] -> ranks [2.0, 1.0, 3.0]
        // Row 1: [0.0, 2.0, 0.0] -> ranks [1.5, 3.0, 1.5]
        assert_eq!(to_f32_vec(&result[0]), vec![2.0, 1.0, 3.0]);
        assert_eq!(to_f32_vec(&result[1]), vec![1.5, 3.0, 1.5]);
    }

    #[test]
    fn test_transposed_output() {
        // Matrix:
        // [1.0, 0.0, 3.0]
        // [0.0, 2.0, 0.0]
        let row_ptr = vec![0, 2, 3];
        let col_indices = vec![0, 2, 1];
        let data = f16_vec(&[1.0, 3.0, 2.0]);

        let result = fast_csr_ranking(&row_ptr, &col_indices, &data, 2, 3, true);

        // Row rankings: [[2.0, 1.0, 3.0], [1.5, 3.0, 1.5]]
        // Transposed structure: outer=columns, inner=rows
        // Col 0: [2.0, 1.5] (row0_col0, row1_col0)
        // Col 1: [1.0, 3.0] (row0_col1, row1_col1)
        // Col 2: [3.0, 1.5] (row0_col2, row1_col2)
        assert_eq!(to_f32_vec(&result[0]), vec![2.0, 1.5]);
        assert_eq!(to_f32_vec(&result[1]), vec![1.0, 3.0]);
        assert_eq!(to_f32_vec(&result[2]), vec![3.0, 1.5]);
    }

    #[test]
    fn test_all_zeros_row() {
        // Matrix
        // [0.0, 0.0, 0.0]
        // [1.0, 2.0, 3.0]
        let row_ptr = vec![0, 0, 3];
        let col_indices = vec![0, 1, 2];
        let data = f16_vec(&[1.0, 2.0, 3.0]);

        let result = fast_csr_ranking(&row_ptr, &col_indices, &data, 2, 3, false);

        // Row 0: all zeros -> all ranks = (1+3)/2 = 2.0
        // Row 1: [1.0, 2.0, 3.0] -> ranks [1.0, 2.0, 3.0]
        assert_eq!(to_f32_vec(&result[0]), vec![2.0, 2.0, 2.0]);
        assert_eq!(to_f32_vec(&result[1]), vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_all_nonzeros_row() {
        // Matrix:
        // [[1.0, 2.0, 3.0]]
        let row_ptr = vec![0, 3];
        let col_indices = vec![0, 1, 2];
        let data = f16_vec(&[1.0, 2.0, 3.0]);

        let result = fast_csr_ranking(&row_ptr, &col_indices, &data, 1, 3, false);

        // No zeros, direct ranking
        assert_eq!(to_f32_vec(&result[0]), vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_tied_values() {
        // Matrix:
        // [[2.0, 0.0, 2.0, 1.0]]
        let row_ptr = vec![0, 3];
        let col_indices = vec![0, 2, 3];
        let data = f16_vec(&[2.0, 2.0, 1.0]);

        let result = fast_csr_ranking(&row_ptr, &col_indices, &data, 1, 4, false);

        // Values: [2.0, 0.0, 2.0, 1.0]
        // Sorted: 0.0(rank1), 1.0(rank2), 2.0(rank3.5), 2.0(rank3.5)
        // Expected: [3.5, 1.0, 3.5, 2.0]
        let expected = [3.5, 1.0, 3.5, 2.0];
        let actual = to_f32_vec(&result[0]);

        for (a, e) in actual.iter().zip(expected.iter()) {
            assert!((a - e).abs() < 0.01, "Expected {}, got {}", e, a);
        }
    }

    #[test]
    fn test_multiple_tied_zeros() {
        // Matrix: [[1.0, 0.0, 0.0, 2.0]]
        let row_ptr = vec![0, 2];
        let col_indices = vec![0, 3];
        let data = f16_vec(&[1.0, 2.0]);

        let result = fast_csr_ranking(&row_ptr, &col_indices, &data, 1, 4, false);

        // Values: [1.0, 0.0, 0.0, 2.0]
        // Two zeros get average rank (1+2)/2 = 1.5
        // Non-zeros: 1.0 gets rank 3, 2.0 gets rank 4
        // Expected: [3.0, 1.5, 1.5, 4.0]
        let expected = [3.0, 1.5, 1.5, 4.0];
        let actual = to_f32_vec(&result[0]);

        for (a, e) in actual.iter().zip(expected.iter()) {
            assert!((a - e).abs() < 0.01, "Expected {}, got {}", e, a);
        }
    }

    #[test]
    fn test_transpose_consistency() {
        // Matrix
        // [1.0, 2.0, 5.0]
        // [3.0, 0.0, 0.0]
        let row_ptr = vec![0, 3, 4];
        let col_indices = vec![0, 1, 2, 0];
        let data = f16_vec(&[1.0, 2.0, 5.0, 3.0]);

        let row_result = fast_csr_ranking(&row_ptr, &col_indices, &data, 2, 3, false);
        let col_result = fast_csr_ranking(&row_ptr, &col_indices, &data, 2, 3, true);

        // Row result:
        // Row 0: [1.0, 2.0, 5.0] -> ranks [1.0, 2.0, 3.0]
        // Row 1: [3.0, 0.0, 0.0] -> ranks [3.0, 1.5, 1.5]

        // Col result should be transpose of row result:
        // Col 0: [1.0, 3.0] (row0_col0, row1_col0)
        // Col 1: [2.0, 1.5] (row0_col1, row1_col1)
        // Col 2: [3.0, 1.5] (row0_col2, row1_col2)

        assert_eq!(row_result.len(), 2); // 2 rows
        assert_eq!(col_result.len(), 3); // 3 cols
        assert_eq!(col_result[0].len(), 2); // 2 rows each

        // Verify transpose relationship
        #[allow(clippy::needless_range_loop)]
        for col_idx in 0..3 {
            for row_idx in 0..2 {
                assert_eq!(
                    row_result[row_idx][col_idx], col_result[col_idx][row_idx],
                    "Mismatch at position ({}, {})",
                    row_idx, col_idx
                );
            }
        }

        // Verify specific expected values
        assert_eq!(to_f32_vec(&row_result[0]), vec![1.0, 2.0, 3.0]);
        assert_eq!(to_f32_vec(&row_result[1]), vec![3.0, 1.5, 1.5]);
        assert_eq!(to_f32_vec(&col_result[0]), vec![1.0, 3.0]);
        assert_eq!(to_f32_vec(&col_result[1]), vec![2.0, 1.5]);
        assert_eq!(to_f32_vec(&col_result[2]), vec![3.0, 1.5]);
    }
}
