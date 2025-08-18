use rand::prelude::*;
use rand_distr::weighted::WeightedIndex;

/// A type alias representing effect size results
///
/// ### Fields
///
/// * `0` - The data in the CSR format
/// * `1` - The column indices
/// * `2` - The row pointers
#[allow(dead_code)]
pub type CscData<T> = (Vec<T>, Vec<usize>, Vec<usize>);

/// Create weighted sparse data resembling single cell counts
///
/// ### Params
///
/// * `n_genes` - Total no of genes
/// * `n_cells` - Total no of cells
/// * `no_genes_exp` - Tuple representing the minimum number and the maximum
///   number of genes expressed per cell
/// * `max_exp` - Maximum expression a given gene can reach. Expression values
///   will be between `1..max_exp`
/// * `seed` - Seed for reproducibility purposes
///
/// ### Returns
///
/// The `CscData` type with the synthetic data.
pub fn create_sparse_sc_data(
    n_genes: usize,
    n_cells: usize,
    no_genes_exp: (usize, usize),
    max_exp: i32,
    seed: usize,
) -> CscData<i32> {
    let weights: Vec<f64> = (1..=n_genes).map(|i| 1.0 / i as f64).collect();
    let total_entries = n_genes * n_cells;

    let mut indptr = vec![0; n_cells + 1];
    let mut indices = Vec::with_capacity(total_entries);
    let mut data = Vec::with_capacity(total_entries);

    for cell_idx in 0..n_cells {
        let mut rng = StdRng::seed_from_u64(seed as u64 + cell_idx as u64);
        let no_cells_expressed = rng.random_range(no_genes_exp.0..=no_genes_exp.1);
        let dist = WeightedIndex::new(&weights).unwrap();

        let start_idx = indices.len();

        for _ in 0..no_cells_expressed {
            let gene_idx = dist.sample(&mut rng);
            let count = rng.random_range(1..=max_exp);

            indices.push(gene_idx);
            data.push(count);
        }

        let end_idx = indices.len();
        let mut gene_count_pairs: Vec<_> = indices[start_idx..end_idx]
            .iter()
            .zip(&data[start_idx..end_idx])
            .map(|(&gene_idx, &count)| (gene_idx, count))
            .collect();

        gene_count_pairs.sort_by_key(|(gene_idx, _)| *gene_idx);

        for (local_idx, (gene_idx, count)) in gene_count_pairs.iter().enumerate() {
            indices[start_idx + local_idx] = *gene_idx;
            data[start_idx + local_idx] = *count;
        }

        indptr[cell_idx + 1] = indices.len();
    }

    (data, indptr, indices)
}
