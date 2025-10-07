use rayon::prelude::*;

use crate::core::data::sparse_io::CsrCellChunk;
use crate::utils::general::flatten_vector;
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
fn rank_f16(vec: &[F16]) -> Vec<f32> {
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

    let mut ranks: Vec<f32> = vec![0.0; n];
    let mut i = 0;
    while i < n {
        let current_value = indexed_values[i].0;
        let start = i;
        while i < n && indexed_values[i].0 == current_value {
            i += 1;
        }
        let avg_rank = (start + i + 1) as f32 / 2.0;
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
/// * `rank_within_rows` - This boolean controls if the ranking happens within
///   cells (for example for AUCell) or across genes (for example for DGE).
///
/// ### Return
///
/// A `Vec<Vec<f32>>` that pending the rank_within_rows represents the ranks
/// across genes or across cells.
pub fn fast_csr_ranking(
    row_ptr: &[usize],
    col_indices: &[u16],
    data: &[F16],
    nrow: usize,
    ncol: usize,
    rank_within_rows: bool,
) -> Vec<Vec<f32>> {
    if rank_within_rows {
        // Rank genes within each cell
        // This is what we are interested in for AUCell type approaches
        (0..nrow)
            .into_par_iter()
            .map(|row_idx| {
                let start = row_ptr[row_idx];
                let end = row_ptr[row_idx + 1];
                let num_nonzeros = end - start;
                let num_zeros = ncol - num_nonzeros;

                if num_nonzeros == 0 {
                    let zero_rank = (1.0 + ncol as f32) / 2.0;
                    return vec![zero_rank; ncol];
                }

                if num_zeros == 0 {
                    let row_data = &data[start..end];
                    return rank_f16(row_data);
                }

                let row_data = &data[start..end];
                let row_cols = &col_indices[start..end];
                let nonzero_ranks = rank_f16(row_data);
                let zero_rank = (1.0 + num_zeros as f32) / 2.0;
                let mut result = vec![zero_rank; ncol];

                for (i, &col) in row_cols.iter().enumerate() {
                    result[col as usize] = nonzero_ranks[i] + num_zeros as f32;
                }

                result
            })
            .collect()
    } else {
        // Rank cells within each gene - build gene-to-cells mapping first
        let mut gene_data: Vec<Vec<(F16, usize)>> = vec![Vec::new(); ncol];

        // Single pass: collect all data per gene
        for row_idx in 0..nrow {
            let start = row_ptr[row_idx];
            let end = row_ptr[row_idx + 1];

            for i in 0..(end - start) {
                let col_idx = col_indices[start + i] as usize;
                gene_data[col_idx].push((data[start + i], row_idx));
            }
        }

        // Rank each gene in parallel
        gene_data
            .into_par_iter()
            .map(|mut values| {
                let num_nonzeros = values.len();
                let num_zeros = nrow - num_nonzeros;

                if num_nonzeros == 0 {
                    let zero_rank = (1.0 + nrow as f32) / 2.0;
                    return vec![zero_rank; nrow];
                }

                values.sort_unstable_by(|a, b| {
                    a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal)
                });

                let zero_rank = (1.0 + num_zeros as f32) / 2.0;
                let mut result = vec![zero_rank; nrow];

                let mut i = 0;
                while i < num_nonzeros {
                    let start_idx = i;
                    let current_value = values[i].0;
                    while i < num_nonzeros && values[i].0 == current_value {
                        i += 1;
                    }
                    let avg_rank = (start_idx + i + 1 + 2 * num_zeros) as f32 / 2.0;
                    for j in start_idx..i {
                        result[values[j].1] = avg_rank;
                    }
                }

                result
            })
            .collect()
    }
}

/// Helper function to rank all cells within a given chunk vector
///
/// ### Params
///
/// * `chunk_vec` - Vector of `CsrCellChunk` to rank.
/// * `no_genes` - Number of represented genes in this data.
/// * `rank_within_rows` - This boolean controls if the ranking happens within
///   cells (for example for AUCell) or across genes (for example for DGE).
///
/// ### Returns
///
/// A `Vec<Vec<f32>>` that pending the rank_within_rows represents the ranks
/// across genes or across cells.
pub fn rank_csr_chunk_vec(
    chunk_vec: Vec<CsrCellChunk>,
    no_genes: usize,
    rank_within_rows: bool,
) -> Vec<Vec<f32>> {
    let no_cells = chunk_vec.len();
    let mut all_data: Vec<Vec<F16>> = Vec::with_capacity(chunk_vec.len());
    let mut all_indices: Vec<Vec<u16>> = Vec::with_capacity(chunk_vec.len());
    let mut indptr: Vec<usize> = Vec::with_capacity(chunk_vec.len() + 1);
    let mut current_indptr = 0_usize;

    indptr.push(current_indptr);

    for chunk in chunk_vec {
        let data_len = chunk.data_norm.len();
        all_data.push(chunk.data_norm);
        all_indices.push(chunk.indices);
        current_indptr += data_len;
        indptr.push(current_indptr);
    }

    let all_data = flatten_vector(all_data);
    let all_indices = flatten_vector(all_indices);

    fast_csr_ranking(
        &indptr,
        &all_indices,
        &all_data,
        no_cells,
        no_genes,
        rank_within_rows,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::traits::F16;

    // Helper to create F16 from f32
    fn f16_vec(values: &[f32]) -> Vec<F16> {
        values.iter().map(|&v| F16::from_f32(v)).collect()
    }

    #[test]
    fn test_simple_row_ranking() {
        // Matrix
        // [1.0, 0.0, 3.0]
        // [0.0, 2.0, 0.0]
        let row_ptr = vec![0, 2, 3];
        let col_indices = vec![0, 2, 1];
        let data = f16_vec(&[1.0, 3.0, 2.0]);

        let result = fast_csr_ranking(&row_ptr, &col_indices, &data, 2, 3, true);

        // Row 0: [1.0, 0.0, 3.0] -> ranks [2.0, 1.0, 3.0]
        // Row 1: [0.0, 2.0, 0.0] -> ranks [1.5, 3.0, 1.5]

        assert_eq!(result[0], vec![2.0, 1.0, 3.0]);
        assert_eq!(result[1], vec![1.5, 3.0, 1.5]);
    }

    #[test]
    fn test_column_ranking() {
        // Matrix:
        // [1.0, 0.0, 3.0]
        // [0.0, 2.0, 0.0]
        let row_ptr = vec![0, 2, 3];
        let col_indices = vec![0, 2, 1];
        let data = f16_vec(&[1.0, 3.0, 2.0]);

        let result = fast_csr_ranking(&row_ptr, &col_indices, &data, 2, 3, false);

        // Now ranking cells within genes:
        // Gene 0: [1.0, 0.0] -> ranks [2.0, 1.0]
        // Gene 1: [0.0, 2.0] -> ranks [1.0, 2.0]
        // Gene 2: [3.0, 0.0] -> ranks [2.0, 1.0]
        assert_eq!(result[0], vec![2.0, 1.0]);
        assert_eq!(result[1], vec![1.0, 2.0]);
        assert_eq!(result[2], vec![2.0, 1.0]);
    }

    #[test]
    fn test_all_zeros_row() {
        // Matrix
        // [0.0, 0.0, 0.0]
        // [1.0, 2.0, 3.0]
        let row_ptr = vec![0, 0, 3];
        let col_indices = vec![0, 1, 2];
        let data = f16_vec(&[1.0, 2.0, 3.0]);

        let result = fast_csr_ranking(&row_ptr, &col_indices, &data, 2, 3, true);

        // Row 0: all zeros -> all ranks = (1+3)/2 = 2.0
        // Row 1: [1.0, 2.0, 3.0] -> ranks [1.0, 2.0, 3.0]
        assert_eq!(result[0], vec![2.0, 2.0, 2.0]);
        assert_eq!(result[1], vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_all_nonzeros_row() {
        // Matrix:
        // [[1.0, 2.0, 3.0]]
        let row_ptr = vec![0, 3];
        let col_indices = vec![0, 1, 2];
        let data = f16_vec(&[1.0, 2.0, 3.0]);

        let result = fast_csr_ranking(&row_ptr, &col_indices, &data, 1, 3, true);

        // No zeros, direct ranking
        assert_eq!(result[0], vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn test_tied_values() {
        // Matrix:
        // [[2.0, 0.0, 2.0, 1.0]]
        let row_ptr = vec![0, 3];
        let col_indices = vec![0, 2, 3];
        let data = f16_vec(&[2.0, 2.0, 1.0]);

        let result = fast_csr_ranking(&row_ptr, &col_indices, &data, 1, 4, true);

        // Values: [2.0, 0.0, 2.0, 1.0]
        // Sorted: 0.0(rank1), 1.0(rank2), 2.0(rank3.5), 2.0(rank3.5)
        // Expected: [3.5, 1.0, 3.5, 2.0]
        let expected = [3.5, 1.0, 3.5, 2.0];
        let actual = &result[0];

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

        let result = fast_csr_ranking(&row_ptr, &col_indices, &data, 1, 4, true);

        // Values: [1.0, 0.0, 0.0, 2.0]
        // Two zeros get average rank (1+2)/2 = 1.5
        // Non-zeros: 1.0 gets rank 3, 2.0 gets rank 4
        // Expected: [3.0, 1.5, 1.5, 4.0]
        let expected = [3.0, 1.5, 1.5, 4.0];
        let actual = &result[0];

        for (a, e) in actual.iter().zip(expected.iter()) {
            assert!((a - e).abs() < 0.01, "Expected {}, got {}", e, a);
        }
    }

    #[test]
    fn test_row_vs_column_ranking() {
        // Matrix
        // [1.0, 2.0, 5.0]
        // [3.0, 0.0, 0.0]
        let row_ptr = vec![0, 3, 4];
        let col_indices = vec![0, 1, 2, 0];
        let data = f16_vec(&[1.0, 2.0, 5.0, 3.0]);

        // Row ranking (rank genes within cells)
        let row_result = fast_csr_ranking(&row_ptr, &col_indices, &data, 2, 3, true);

        // Column ranking (rank cells within genes)
        let col_result = fast_csr_ranking(&row_ptr, &col_indices, &data, 2, 3, false);

        // Row result (rank genes within each cell):
        // Row 0: [1.0, 2.0, 5.0] -> ranks [1.0, 2.0, 3.0]
        // Row 1: [3.0, 0.0, 0.0] -> ranks [3.0, 1.5, 1.5]
        assert_eq!(row_result.len(), 2);
        assert_eq!(row_result[0], vec![1.0, 2.0, 3.0]);
        assert_eq!(row_result[1], vec![3.0, 1.5, 1.5]);

        // Column result (rank cells within each gene):
        // Gene 0: [1.0, 3.0] -> ranks [1.0, 2.0]
        // Gene 1: [2.0, 0.0] -> ranks [2.0, 1.0]
        // Gene 2: [5.0, 0.0] -> ranks [2.0, 1.0]
        assert_eq!(col_result.len(), 3);
        assert_eq!(col_result[0], vec![1.0, 2.0]);
        assert_eq!(col_result[1], vec![2.0, 1.0]);
        assert_eq!(col_result[2], vec![2.0, 1.0]);
    }

    #[test]
    fn test_column_ranking_with_ties() {
        // Matrix
        // [2.0, 1.0]
        // [0.0, 1.0]
        // [2.0, 0.0]
        let row_ptr = vec![0, 2, 3, 4];
        let col_indices = vec![0, 1, 1, 0];
        let data = f16_vec(&[2.0, 1.0, 1.0, 2.0]);

        let result = fast_csr_ranking(&row_ptr, &col_indices, &data, 3, 2, false);

        // Gene 0: [2.0, 0.0, 2.0] -> tied 2.0s at positions 0 and 2
        //   Sorted: 0.0(rank1), 2.0(rank2.5), 2.0(rank2.5)
        //   Result: [2.5, 1.0, 2.5]
        // Gene 1: [1.0, 1.0, 0.0] -> tied 1.0s at positions 0 and 1
        //   Sorted: 0.0(rank1), 1.0(rank2.5), 1.0(rank2.5)
        //   Result: [2.5, 2.5, 1.0]

        let gene0_actual = &result[0];
        let gene1_actual = &result[1];

        assert!((gene0_actual[0] - 2.5).abs() < 0.01);
        assert!((gene0_actual[1] - 1.0).abs() < 0.01);
        assert!((gene0_actual[2] - 2.5).abs() < 0.01);

        assert!((gene1_actual[0] - 2.5).abs() < 0.01);
        assert!((gene1_actual[1] - 2.5).abs() < 0.01);
        assert!((gene1_actual[2] - 1.0).abs() < 0.01);
    }
}
