use extendr_api::prelude::*;

use crate::helpers::synthetic_data::*;
use crate::utils::r_rust_interface::{faer_to_r_matrix, r_matrix_to_faer};

/// @export
#[extendr]
fn rs_generate_bulk_rnaseq(
    num_samples: usize,
    num_genes: usize,
    seed: usize,
    add_modules: bool,
    module_sizes: Option<Vec<i32>>,
) -> List {
    // r cannot deal with usize
    let module_sizes: Vec<usize> = module_sizes
        .unwrap_or(vec![300, 250, 200, 300, 500])
        .iter()
        .map(|x| *x as usize)
        .collect();

    let data = generate_bulk_rnaseq(
        num_samples,
        num_genes,
        seed as u64,
        add_modules,
        Some(module_sizes),
    );

    let matrix = faer_to_r_matrix(data.count_matrix.as_ref());
    let module_membership: Vec<i32> = data.gene_modules.iter().map(|x| *x as i32).collect();

    list!(counts = matrix, module_membership = module_membership)
}

/// @export
#[extendr]
fn rs_simulate_dropouts(
    count_mat: RMatrix<f64>,
    dropout_midpoint: f64,
    dropout_shape: f64,
    global_sparsity: f64,
    seed: usize,
) -> extendr_api::RArray<f64, [usize; 2]> {
    let data = r_matrix_to_faer(&count_mat);

    let sparse_data = simulate_dropouts(
        &data,
        dropout_midpoint,
        dropout_shape,
        global_sparsity,
        seed as u64,
    );

    faer_to_r_matrix(sparse_data.as_ref())
}

extendr_module! {
    mod r_synthetic_data;
    fn rs_generate_bulk_rnaseq;
    fn rs_simulate_dropouts;
}
