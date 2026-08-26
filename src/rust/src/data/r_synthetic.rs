use extendr_api::prelude::*;

use bixverse_rs::core::synthetic_data::*;
use bixverse_rs::prelude::*;
use bixverse_rs::single_cell::sc_data::sc_synthetic_data::*;

/////////////
// extendR //
/////////////

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
    // -- dialogue --
    fn rs_synthetic_sc_dialogue_data;
}

//////////
// Bulk //
//////////

/// Generation of bulkRNAseq-like data with optional correlation structure
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Function generates synthetic bulkRNAseq data with heteroskedasticity (lowly
/// expressed genes show higher variance) and optional co-expression modules
/// planted on a latent factor. Alongside the counts it returns the ground truth
/// (module membership, hub genes, per-gene loadings and the latent factors), so
/// downstream methods can be scored against what was actually simulated.
///
/// @param synthetic_params List. The synthetic data parameters, see
/// [bixverse::params_synthetic_bulk_rnaseq()]. Expected elements are
/// `num_samples`, `num_genes`, `module_sizes` (integer vector, empty means no
/// modules), `generator` (one of `c("hub_modular", "modular",
/// "non_negative_factor", "non_gaussian_factor")`), `seed`,
/// `mean_exp_gamma_shape`, `mean_exp_gamma_scale`, `disp_intercept`,
/// `disp_slope`, `noise_std`, `factor_std`, `factor_shape`, `factor_scale`,
/// `loading_mu`, `loading_sigma` and `hub_percentile`.
///
/// @return List with the following elements
/// \itemize{
///     \item counts The matrix of simulated counts. Rows are genes, columns
///     are samples.
///     \item module_membership Vector defining the module membership. `0` is
///     background, `1..K` the module identifier.
///     \item module_hubs 1-indexed positions of the genes flagged as hubs.
///     Empty for the `"modular"` generator, which plants no hubs.
///     \item loadings Per-gene loading on its module's latent factor. `0` for
///     background genes.
///     \item module_factors The latent factor matrix. Rows are modules,
///     columns are samples.
/// }
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_generate_bulk_rnaseq(synthetic_params: List) -> extendr_api::Result<List> {
    let params: SyntheticRnaSeqParams<f64> = SyntheticRnaSeqParams::from_r_list(synthetic_params)?;

    let data: SyntheticRnaSeqData<f64> = generate_bulk_rnaseq(&params).to_extendr()?;

    Ok(list!(
        counts = faer_to_r_matrix(data.count_matrix.as_ref()),
        module_membership = data.gene_modules.r_int_convert(),
        // Rust hands back 0-indexed gene positions; R wants 1-indexed
        module_hubs = data.module_hubs.r_int_convert_shift(),
        loadings = data.loadings,
        module_factors = faer_to_r_matrix(data.module_factors.as_ref())
    ))
}

/// Sparsify bulkRNAseq like data
///
/// @description
/// `r lifecycle::badge("experimental")`
/// This function takes in a (raw) count matrix (for example from the synthetic
/// data in bixverse) and applies Splatter-style sequencing-depth dropout to it.
/// Per sample a size factor `s_j ~ LogNormal(0, capture_efficiency_sigma)` is
/// drawn, giving a target library size of `target_library_size * s_j`. Each
/// gene is then binomially thinned to approach that target. Retention
/// probability is capped at 1, so samples below their target are left alone
/// rather than upsampled.
///
/// @param count_mat Numerical matrix. Original numeric matrix. Rows are genes,
/// columns are samples.
/// @param sparsity_params List. The sparsity parameters, see
/// [bixverse::params_bulk_sparsity()]. Expected elements are `strategy`,
/// `target_library_size`, `capture_efficiency_sigma` and `seed`.
///
/// @return The sparsified matrix based on the provided parameters.
///
/// @export
///
/// @references Zappia, et al., Genome Biol, 2017
///
/// @keywords internal
#[extendr]
fn rs_simulate_dropouts(
    count_mat: RMatrix<f64>,
    sparsity_params: List,
) -> extendr_api::Result<RArray<f64, 2>> {
    let data = r_matrix_to_faer(&count_mat);

    let params: SparsityParams<f64> = SparsityParams::from_r_list(sparsity_params)?;

    let sparse_data = apply_dropout(data.as_ref(), &params).to_extendr()?;

    Ok(faer_to_r_matrix(sparse_data.as_ref()))
}

/////////////////
// Single cell //
/////////////////

/// Generates synthetic data for single cell
///
/// @description
/// `r lifecycle::badge("experimental")`
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
///
/// @keywords internal
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
/// `r lifecycle::badge("experimental")`
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
///
/// @keywords internal
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
/// `r lifecycle::badge("experimental")`
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
///
/// @export
///
/// @keywords internal
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

//////////////
// DIALOGUE //
//////////////

/// Generates synthetic data with a planted multicellular programme
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Builds the fixture DIALOGUE is tested against. Every cell type gets its own
/// noise and its own sample-level nuisance factors; only the first feature
/// column and the planted genes carry the shared per-sample latent, so anything
/// found beyond that is spurious. Cells are laid out contiguously by cell type
/// and, within a cell type, by sample.
///
/// The counts are a scaled copy of the normalised layer rather than a draw from
/// a count model. That keeps the planted signal clean, which is the point of a
/// fixture.
///
/// @param n_samples Integer. Samples the experiment spans. DIALOGUE needs at
/// least 5.
/// @param cells_per_sample Integer. Cells per sample per cell type.
/// @param n_cell_types Integer. Number of cell types. Must be at least 2.
/// @param n_features Integer. Feature columns per cell type. Must be at least
/// 2.
/// @param n_sample_features Integer. Feature columns carrying a per-sample
/// component. The first of those is the shared programme, the rest are
/// cell-type-specific nuisance; anything past this count is pure noise.
/// @param n_genes Integer. Number of genes.
/// @param n_planted Integer. Planted genes per cell type. Cell type `t` owns
/// genes `t * n_planted` to `(t + 1) * n_planted - 1` (0-indexed), so the
/// blocks have to fit into `n_genes`.
/// @param seed Integer. Random seed for reproducibility.
///
/// @return A list with the following items.
/// \itemize{
///   \item data - The synthetic raw counts, CSR over cells.
///   \item indptr - The index pointers of the cells.
///   \item indices - The gene indices for the given cells.
///   \item nrow - Number of cells.
///   \item ncol - Number of genes.
///   \item cell_type_indices - List of integer vectors. 0-indexed(!) global
///   cell positions per cell type.
///   \item features - List of numeric matrices, one per cell type, rows
///   aligned to `cell_type_indices`.
///   \item sample_ids - Integer vector. 0-indexed(!) sample code per cell.
///   \item quality - Numeric vector. Quality covariate per cell. Pure noise.
///   \item latent - Numeric vector. The per-sample latent the planted
///   programme follows.
///   \item planted - List of integer vectors. 0-indexed(!) planted gene
///   positions per cell type.
/// }
///
/// @references Jerby-Arnon & Regev, Nature Biotechnology, 2022
///
/// @export
///
/// @keywords internal
#[allow(clippy::too_many_arguments)]
#[extendr]
fn rs_synthetic_sc_dialogue_data(
    n_samples: usize,
    cells_per_sample: usize,
    n_cell_types: usize,
    n_features: usize,
    n_sample_features: usize,
    n_genes: usize,
    n_planted: usize,
    seed: usize,
) -> extendr_api::Result<List> {
    let params = DialogueSyntheticParams::new(
        n_samples,
        cells_per_sample,
        n_cell_types,
        n_features,
        n_sample_features,
        n_genes,
        n_planted,
    );

    let synthetic = create_dialogue_synthetic_data(&params, seed as u64).to_extendr()?;

    // The generator works gene-major; R wants a `dgRMatrix`, so hand back the
    // cell-major view of the same matrix.
    let csr = synthetic.matrix.transform();

    let features = List::from_values(
        synthetic
            .features
            .iter()
            .map(|mat| faer_to_r_matrix(mat.as_ref())),
    );
    let cell_type_indices = List::from_values(
        synthetic
            .cell_type_indices
            .into_iter()
            .map(|cells| cells.r_int_convert()),
    );
    let planted = List::from_values(
        synthetic
            .planted
            .into_iter()
            .map(|genes| genes.r_int_convert()),
    );

    Ok(list!(
        data = csr.data,
        indptr = csr.indptr,
        indices = csr.indices,
        nrow = csr.shape.0,
        ncol = csr.shape.1,
        cell_type_indices = cell_type_indices,
        features = features,
        sample_ids = synthetic.sample_ids.r_int_convert(),
        quality = synthetic.quality,
        latent = synthetic.latent,
        planted = planted
    ))
}
