use extendr_api::prelude::*;

use crate::helpers::structs_sparse::{CscData, CsrData};
use crate::helpers::synthetic_data::*;

/// Generate synthetic single cell data (Seurat type)
///
/// @description This function generates pseudo data to test single cell
/// functions in form of the Seurat version, with cells = columns and genes =
/// rows. The data is CSC type.
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
/// @return The list with the synthetic data with the following items:
///  \itemize{
///   \item data - The synthetic counts
///   \item col_ptrs - The column pointers
///   \item row_indices - The row indices
/// }
///
/// @export
#[extendr]
fn rs_synthetic_sc_data_csc(
    n_genes: usize,
    n_cells: usize,
    min_genes: usize,
    max_genes: usize,
    max_exp: i32,
    seed: usize,
) -> List {
    let synthetic_data: CscData<i32> =
        create_sparse_csc_data(n_genes, n_cells, (min_genes, max_genes), max_exp, seed);

    list!(
        data = synthetic_data.0,
        col_ptrs = synthetic_data.1,
        row_indices = synthetic_data.2
    )
}

/// Generate synthetic single cell data (h5ad type)
///
/// @description This function generates pseudo data to test single cell
/// functions in form of the h5ad version, with cells = rows and genes =
/// columns. The data is CSR type.
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
/// @return The list with the synthetic data with the following items:
///  \itemize{
///   \item data - The synthetic counts
///   \item row_ptrs - The row pointers
///   \item col_indices - The column indices
/// }
///
/// @export
#[extendr]
fn rs_synthetic_sc_data_csr(
    n_genes: usize,
    n_cells: usize,
    min_genes: usize,
    max_genes: usize,
    max_exp: i32,
    seed: usize,
) -> List {
    let synthetic_data: CsrData<i32> =
        create_sparse_csr_data(n_genes, n_cells, (min_genes, max_genes), max_exp, seed);

    list!(
        data = synthetic_data.0,
        row_ptrs = synthetic_data.1,
        col_indices = synthetic_data.2
    )
}

extendr_module! {
    mod r_synthetic_data;
    fn rs_synthetic_sc_data_csc;
    fn rs_synthetic_sc_data_csr;
}
