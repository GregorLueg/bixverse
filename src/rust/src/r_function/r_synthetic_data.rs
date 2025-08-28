use extendr_api::prelude::*;

use crate::helpers::synthetic_data::*;
use crate::utils::r_rust_interface::{faer_to_r_matrix, r_matrix_to_faer};

/// Generation of bulkRNAseq-like data with optional correlation structure
///
/// @description
/// Function generates synthetic bulkRNAseq data with heteroskedasticity (lowly
/// expressed genes show higher variance) and can optionally add correlation
/// structures for testing purposes.
///
/// @param num_samples Integer. Number of samples to simulate.
/// @param num_genes Integer. Number of genes to simulate.
/// @param seed. Integer. Seed for reproducibility.
/// @param add_modules Boolean. Shall correlation structures be added to the
/// data.
/// @param module_size `NULL` or vector of sizes of the gene modules. When
/// `NULL` defaults to `c(300, 250, 200, 300, 500)`. Warning! The sum of this
/// vector must be ≤ num_genes!
///
/// @return List with the following elements
/// \itemize{
///     \item counts The matrix of simulated counts.
///     \item module_membership Vector defining the module membership.
/// }
///
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

/// Sparsify bulkRNAseq like data
///
/// @description
/// This function takes in a (raw) count matrix (for example from the synthetic
/// data in bixverse) and applies sparsification to it based on two possible
/// functions:
///
/// **Logistic function:**
///
/// With dropout probability defined as:
///
/// ```
/// P(dropout) = clamp(1 / (1 + exp(shape * (ln(exp+1) - ln(midpoint+1)))), 0.3, 0.8) * (1 - global_sparsity) + global_sparsity
/// ```
///
/// with the following characteristics:
///
/// - Plateaus at global_sparsity dropout for high expression genes
/// - Partial dropout preserves count structure via binomial thinning
/// - Good for preserving variance-mean relationships
///
/// **Power Decay function:**
///
/// With dropout probability defined as:
///
/// ```
/// P(dropout) = (midpoint / (exp + midpoint))^power * scale_factor * (1 - global_sparsity) + global_sparsity
/// ```
///
/// with the following characteristics:
///
/// - No plateau - high expression genes get substantial dropout
/// - Complete dropout only (no partial dropout)
/// - More uniform dropout across expression range
///
/// @param count_mat Numerical matrix. Original numeric matrix.
/// @param dropout_function String. One of `c("log", "powerdecay")`. Defines
/// which function will be used to induce the sparsity.
/// @param dropout_midpoint Numeric. Controls the midpoint parameter of the
/// logistic and power decay function.
/// @param dropout_shape Numeric. Controls the shape parameter of the logistic
/// function.
/// @param power_factor Numeric. Controls the power factor of the power decay
/// function.
/// @param global_sparsity Numeric. The global sparsity parameter.
/// @param seed Integer. Seed for reproducibility.
///
/// @return The sparsified matrix based on the provided parameters.
///
/// @export
#[extendr]
fn rs_simulate_dropouts(
    count_mat: RMatrix<f64>,
    dropout_function: String,
    dropout_midpoint: f64,
    dropout_shape: f64,
    power_factor: f64,
    global_sparsity: f64,
    seed: usize,
) -> extendr_api::Result<extendr_api::RArray<f64, [usize; 2]>> {
    let data = r_matrix_to_faer(&count_mat);

    let dropout_type = parse_sparsification(&dropout_function)
        .ok_or_else(|| format!("Invalid dropout_function type: {}", dropout_function))?;

    let sparse_data = match dropout_type {
        SparsityFunction::Logistic => simulate_dropouts_logistic(
            &data,
            dropout_midpoint,
            dropout_shape,
            global_sparsity,
            seed as u64,
        ),
        SparsityFunction::PowerDecay => simulate_dropouts_power_decay(
            &data,
            dropout_midpoint,
            power_factor,
            global_sparsity,
            seed as u64,
        ),
    };

    Ok(faer_to_r_matrix(sparse_data.as_ref()))
}

extendr_module! {
    mod r_synthetic_data;
    fn rs_generate_bulk_rnaseq;
    fn rs_simulate_dropouts;
}
