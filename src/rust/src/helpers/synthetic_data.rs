use rand::prelude::*;
use rand_distr::weighted::WeightedAliasIndex;

use crate::helpers::structs_sparse::{CscData, CsrData};

/// Create weighted sparse data resembling single cell counts in CSC
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
pub fn create_sparse_csc_data(
    n_genes: usize,
    n_cells: usize,
    no_cells_exp: (usize, usize),
    max_exp: i32,
    seed: usize,
) -> CscData<i32> {
    let weights: Vec<f64> = (1..=n_genes).map(|i| 1.0 / i as f64).collect();

    let avg_cells = (no_cells_exp.0 + no_cells_exp.1) / 2;
    let estimated_total = n_genes * avg_cells;

    let mut indptr = Vec::with_capacity(n_genes + 1);
    let mut indices = Vec::with_capacity(estimated_total);
    let mut data = Vec::with_capacity(estimated_total);
    indptr.push(0);

    let mut temp_vec = Vec::with_capacity(no_cells_exp.1);

    #[allow(clippy::needless_range_loop)]
    for gene_idx in 0..n_genes {
        let mut rng = StdRng::seed_from_u64(seed as u64 + gene_idx as u64);

        // Trick thanks to Claude
        let weight_factor = weights[gene_idx];
        let scaled_range = (
            (no_cells_exp.0 as f64 * weight_factor * 10.0) as usize,
            (no_cells_exp.1 as f64 * weight_factor * 10.0).min(n_cells as f64) as usize,
        );
        let no_cells_expressing = rng.random_range(scaled_range.0.max(1)..=scaled_range.1.max(1));

        temp_vec.clear();

        for _ in 0..no_cells_expressing {
            let cell_idx = rng.random_range(0..n_cells);
            let count = rng.random_range(1..=max_exp);
            temp_vec.push((cell_idx, count));
        }

        temp_vec.sort_unstable_by_key(|(cell_idx, _)| *cell_idx);

        // remove duplicates by summing counts
        temp_vec.dedup_by(|a, b| {
            if a.0 == b.0 {
                b.1 += a.1;
                true
            } else {
                false
            }
        });

        for (cell_idx, count) in temp_vec.iter() {
            indices.push(*cell_idx); // Cell indices for CSC
            data.push(*count);
        }

        indptr.push(indices.len());
    }

    (data, indptr, indices, None)
}

/// Create weighted sparse data resembling single cell counts in CSR format
///
/// ### Params
///
/// * `n_genes` - Total no of genes  
/// * `n_cells` - Total no of cells
/// * `no_genes_exp` - Tuple representing the min and max number of genes expressed
///   per cell
/// * `max_exp` - Maximum expression a given gene can reach. Expression values will
///   be between `1..max_exp`
/// * `seed` - Seed for reproducibility purposes
///
/// ### Returns
///
/// The `CsrData` type with the synthetic data (cells as rows, genes as columns).
pub fn create_sparse_csr_data(
    n_genes: usize,
    n_cells: usize,
    no_genes_exp: (usize, usize),
    max_exp: i32,
    seed: usize,
) -> CsrData<i32> {
    let weights: Vec<f64> = (1..=n_genes).map(|i| 1.0 / i as f64).collect();
    let alias = WeightedAliasIndex::new(weights).unwrap();

    let avg_genes = (no_genes_exp.0 + no_genes_exp.1) / 2;
    let estimated_total = n_cells * avg_genes;

    let mut indptr = Vec::with_capacity(n_cells + 1);
    let mut indices = Vec::with_capacity(estimated_total);
    let mut data = Vec::with_capacity(estimated_total);
    indptr.push(0);

    let mut temp_vec = Vec::with_capacity(no_genes_exp.1);

    for cell_idx in 0..n_cells {
        let mut rng = StdRng::seed_from_u64(seed as u64 + cell_idx as u64);
        let no_genes_expressed = rng.random_range(no_genes_exp.0..=no_genes_exp.1);

        temp_vec.clear();

        for _ in 0..no_genes_expressed {
            let gene_idx = alias.sample(&mut rng);
            let count = rng.random_range(1..=max_exp);
            temp_vec.push((gene_idx, count));
        }

        // Sort by gene index
        temp_vec.sort_unstable_by_key(|(gene_idx, _)| *gene_idx);

        for (gene_idx, count) in temp_vec.iter() {
            indices.push(*gene_idx);
            data.push(*count);
        }

        indptr.push(indices.len());
    }

    (data, indptr, indices, None)
}
