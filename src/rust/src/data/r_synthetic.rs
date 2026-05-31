use extendr_api::prelude::*;

use bixverse_rs::core::synthetic_data::*;
use bixverse_rs::prelude::*;
use bixverse_rs::single_cell::sc_data::sc_synthetic_data::*;

extendr_module! {
    mod r_synthetic;
    // bulk
    fn rs_generate_bulk_rnaseq;
    fn rs_simulate_dropouts;
    // single cell
    // -- rna --
    fn rs_synthetic_sc_data_with_cell_types;
    fn rs_sample_ids_for_cell_types;
    // -- adt --
    fn rs_synthetic_sc_adt_with_cell_types;
}

//////////
// Bulk //
//////////

/// Generation of bulkRNAseq-like data with optional correlation structure
///
/// @description
/// Function generates synthetic bulkRNAseq data with heteroskedasticity (lowly
/// expressed genes show higher variance) and can optionally add correlation
/// structures for testing purposes.
///
/// @param num_samples Integer. Number of samples to simulate.
/// @param num_genes Integer. Number of genes to simulate.
/// @param seed Integer. Seed for reproducibility.
/// @param add_modules Boolean. Shall correlation structures be added to the
/// data.
/// @param module_sizes `NULL` or vector of sizes of the gene modules. When
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

    let data: SyntheticRnaSeqData<f64> = generate_bulk_rnaseq(
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
) -> extendr_api::Result<RArray<f64, 2>> {
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

/////////////////
// Single cell //
/////////////////

/// Generates synthetic data for single cell
///
/// @description
/// Helper function to generate synthetic single cell data with optional
/// bathc effects and sample bias.
///
/// @param n_cells Integer. Number of cells to generate.
/// @param n_genes Integer. Number of genes to generate.
/// @param n_batches Integer. Number of the batches to generated.
/// @param n_samples Optional integer. Shall the cells be distributed over
/// `n_samples` samples.
/// @param cell_configs A nested list that indicates which gene indices
/// are markers for which cell.
/// @param batch_effect_strength String. One of `c("strong", "medium", "low")`.
/// Defines the strength of the added batch effect.
/// @param sample_bias Optional string. One of
/// `c("even", "slightly_uneven", "very_uneven")`
/// @param seed Integer. Random seed for reproducibility.
///
/// @return A list with the following items.
/// \itemize{
///   \item data - The synthetic raw counts.
///   \item indptr - The index pointers of the cells.
///   \item indices - The indices of the genes for the given cells.
///   \item nrow - Number of rows.
///   \item ncol - Number of columns
///   \item cell_type_indices - Vector indicating which cell type this is.
///   \item batch_indices - Vector indicating the batch.
///   \item sample_indices - Optional sample indices if asked for.
/// }
///
/// @export
#[allow(clippy::too_many_arguments)]
#[extendr]
fn rs_synthetic_sc_data_with_cell_types(
    n_cells: usize,
    n_genes: usize,
    n_batches: usize,
    n_samples: Option<usize>,
    cell_configs: List,
    batch_effect_strength: String,
    sample_bias: Option<String>,
    seed: usize,
) -> extendr_api::Result<List> {
    let mut cell_configs_vec = Vec::with_capacity(cell_configs.len());
    for i in 0..cell_configs.len() {
        let elem_i = cell_configs.elt(i)?;
        let list_i = elem_i
            .as_list()
            .ok_or_else(|| extendr_api::Error::Other("Expected list".into()))?;
        let cell_config = CellTypeConfig::from_r_list(list_i)?;
        cell_configs_vec.push(cell_config);
    }

    let synthetic_data: (CompressedSparseData2<u32>, Vec<usize>, Vec<usize>) =
        create_celltype_sparse_csr_data(
            n_cells,
            n_genes,
            cell_configs_vec,
            n_batches,
            &batch_effect_strength,
            seed,
        );

    match (n_samples, sample_bias) {
        (Some(n_samp), Some(bias_str)) => {
            let bias = parse_sample_bias(&bias_str)
                .ok_or_else(|| extendr_api::Error::Other("Invalid sample_bias value".into()))?;
            let sample_labels = generate_sample_labels(&synthetic_data.1, n_samp, &bias, seed);
            Ok(list!(
                data = synthetic_data.0.data,
                indptr = synthetic_data.0.indptr,
                indices = synthetic_data.0.indices,
                nrow = synthetic_data.0.shape.0,
                ncol = synthetic_data.0.shape.1,
                cell_type_indices = synthetic_data.1.r_int_convert(),
                batch_indices = synthetic_data.2.r_int_convert(),
                sample_indices = sample_labels.r_int_convert()
            ))
        }
        _ => Ok(list!(
            data = synthetic_data.0.data,
            indptr = synthetic_data.0.indptr,
            indices = synthetic_data.0.indices,
            nrow = synthetic_data.0.shape.0,
            ncol = synthetic_data.0.shape.1,
            cell_type_indices = synthetic_data.1.r_int_convert(),
            batch_indices = synthetic_data.2.r_int_convert(),
            sample_indices = NULL
        )),
    }
}

/// Generates synthetic ADT counts with defined cell types
///
/// @description
/// Generates a dense cells x proteins matrix of synthetic raw ADT counts for
/// testing. Proteins are assigned roles by column index: marker proteins are
/// elevated in their owning cell type and sit at background elsewhere, isotype
/// controls only ever carry background, and any column named as neither is a
/// generic background-only protein. Counts follow a negative-binomial draw with
/// an additive background plus per-cell-type signal, a per-cell capture
/// efficiency factor, and an optional per-batch staining multiplier. Cell type
/// and batch assignment match `rs_synthetic_sc_with_cell_types()` cell-for-cell
/// for matched inputs, so RNA and ADT can be paired for multi-modal tests.
///
/// @param n_cells Integer. Number of cells (matrix rows).
/// @param n_proteins Integer. Total number of proteins (matrix columns). Must
/// be large enough to cover every marker and isotype index supplied.
/// @param n_batches Integer. Number of batches. Batch 0 is unperturbed; further
/// batches receive a per-protein staining multiplier.
/// @param isotype_controls Integer vector. The 0-based column indices that are
/// isotype controls. These are forced to background only, even if they also
/// appear in a cell type's markers.
/// @param cell_configs List. One element per cell type, each a list with a
/// `marker_genes` integer vector of 0-based marker column indices for that
/// cell type.
/// @param batch_effect_strength String. One of `c("weak", "medium", "strong")`.
/// Controls the spread of the per-batch staining multiplier. Unrecognised
/// values fall back to `"strong"`.
/// @param seed Integer. For reproducibility.
///
/// @return A list with the following items.
/// \itemize{
///   \item data - Integer vector. The counts in row-major order, length
///   `n_cells * n_proteins` (cell-major: all proteins of cell 0, then cell 1).
///   \item cell_type_indices - Integer vector of length `n_cells`. The 0-based
///   cell type assigned to each cell.
///   \item batch_indices - Integer vector of length `n_cells`. The 0-based
///   batch assigned to each cell.
/// }
///
/// @export
#[allow(clippy::too_many_arguments)]
#[extendr]
fn rs_synthetic_sc_adt_with_cell_types(
    n_cells: usize,
    n_proteins: usize,
    n_batches: usize,
    isotype_controls: Vec<i32>,
    cell_configs: List,
    batch_effect_strength: String,
    seed: usize,
) -> extendr_api::Result<List> {
    let isotype_controls: Vec<usize> = isotype_controls.r_int_convert();
    let mut cell_configs_vec = Vec::with_capacity(cell_configs.len());
    for i in 0..cell_configs.len() {
        let elem_i = cell_configs.elt(i)?;
        let list_i = elem_i
            .as_list()
            .ok_or_else(|| extendr_api::Error::Other("Expected list".into()))?;
        let cell_config = CellTypeConfig::from_r_list(list_i)?;
        cell_configs_vec.push(cell_config);
    }

    let synthetic_data: (Vec<u32>, Vec<usize>, Vec<usize>) = create_adt_synthetic_data(
        n_cells,
        n_proteins,
        cell_configs_vec,
        isotype_controls,
        n_batches,
        &batch_effect_strength,
        seed,
    );

    Ok(list!(
        data = synthetic_data.0.r_int_convert(),
        cell_type_indices = synthetic_data.1.r_int_convert(),
        batch_indices = synthetic_data.2.r_int_convert()
    ))
}

/// Helper function to generate sample identifiers based on cells
///
/// @description
/// Extract out of `rs_synthetic_sc_data_with_cell_types()` to quickly iterate
/// over different sample to cell type patterns
///
/// @param cell_type_indices Integer vector. Each integer represents a cell
/// type.
/// @param n_samples Integer. Number of different sample ids to generate.
/// @param sample_bias String. One of
/// `c("even", "slightly_uneven", "very_uneven")`. Determins the cell type
/// to sample id associations.
/// @param seed Integer. Random seed for reproducibility.
///
/// @returns An integer vector representing the samples.
#[extendr]
fn rs_sample_ids_for_cell_types(
    cell_type_indices: &[i32],
    n_samples: usize,
    sample_bias: String,
    seed: usize,
) -> extendr_api::Result<Vec<i32>> {
    let cell_type_indices = cell_type_indices.r_int_convert();
    let bias = parse_sample_bias(&sample_bias)
        .ok_or_else(|| extendr_api::Error::Other("Invalid sample_bias value".into()))?;
    let sample_labels = generate_sample_labels(&cell_type_indices, n_samples, &bias, seed);

    Ok(sample_labels.r_int_convert())
}
