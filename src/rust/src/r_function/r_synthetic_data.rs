use extendr_api::prelude::*;

use crate::helpers::synthetic_data::*;

/// Generate synthetic single cell data
///
/// @description This function generates pseudo data to test single cell
/// functions
///
/// @param n_genes Integer. Number of genes you wish to have in the synthetic
/// data.
/// @param n_cells Integer. Number of cells you wish to have in the synthetic
/// data.
/// @param min_genes Integer. Minimum number of genes expressed per cell.
/// @param max_genes Integer. Maximum number of genes expressed per cell.
/// @param max_exp Upper bound in terms of expression. Expression values will be
/// sampled from `1:max_exp`.
/// @param seed Integer. Seed for reproducibility purposes.
///
/// @export
#[extendr]
fn rs_synthetic_sc_data(
    n_genes: usize,
    n_cells: usize,
    min_genes: usize,
    max_genes: usize,
    max_exp: i32,
    seed: usize,
) -> List {
    let synthetic_data =
        create_sparse_sc_data(n_genes, n_cells, (min_genes, max_genes), max_exp, seed);

    list!(
        data = synthetic_data.0,
        col_ptrs = synthetic_data.1,
        row_indices = synthetic_data.2
    )
}

extendr_module! {
    mod r_synthetic_data;
    fn rs_synthetic_sc_data;
}
