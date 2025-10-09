use extendr_api::*;
use faer::Mat;
use indexmap::IndexSet;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use std::time::Instant;
use thousands::Separable;

use crate::core::base::loess::*;
use crate::core::base::pca_svd::{randomised_svd_f32, RandomSvdResults};
use crate::core::data::sparse_io::*;

////////////////
// Structures //
////////////////

/// Structure to store QC information on cells
///
/// ### Fields
///
/// * `cell_indices` - Indices of which cells to keep.
/// * `lib_size` - Optional library size of the cells.
/// * `no_genes` - Optional number of genes of the cells.
#[derive(Clone, Debug)]
#[allow(dead_code)]
pub struct CellQuality {
    pub cell_indices: Vec<usize>,
    pub gene_indices: Vec<usize>,
    pub lib_size: Vec<usize>,
    pub no_genes: Vec<usize>,
}

impl CellQuality {
    /// Update the internal cell indices
    ///
    /// ### Params
    ///
    /// * `cell_indices` - Vector of cell indices to keep
    pub fn set_cell_indices(&mut self, cell_indices: &[usize]) {
        self.cell_indices = cell_indices.to_vec();
    }

    /// Update the internal gene indices
    ///
    /// ### Params
    ///
    /// * `gene_indices` - Vector of gene indices to keep
    pub fn set_gene_indices(&mut self, gene_indices: &[usize]) {
        self.gene_indices = gene_indices.to_vec();
    }
}

/// Structure that stores minimum QC thresholds/info for single cell
///
/// ### Fields
///
/// * `min_unique_genes` - Minimum number of unique genes per cell/spot.
/// * `min_lib_size` - Minimum library size per cell/spot.
/// * `min_cells` - Minimum cells per gene.
/// * `target_size` - Target size for normalisation.
#[derive(Clone, Debug)]
pub struct MinCellQuality {
    pub min_unique_genes: usize,
    pub min_lib_size: usize,
    pub min_cells: usize,
    pub target_size: f32,
}

impl MinCellQuality {
    /// Generate the MinCellQuality params from an R list
    ///
    /// or default to sensible defaults
    ///
    /// ### Params
    ///
    /// * `r_list` - R list with the parameters
    ///
    /// ### Returns
    ///
    /// Self with the specified parameters.
    pub fn from_r_list(r_list: List) -> Self {
        let min_qc = r_list.into_hashmap();

        let min_unique_genes = min_qc
            .get("min_unique_genes")
            .and_then(|v| v.as_integer())
            .unwrap_or(100) as usize;

        let min_lib_size = min_qc
            .get("min_lib_size")
            .and_then(|v| v.as_integer())
            .unwrap_or(250) as usize;

        let min_cells = min_qc
            .get("min_cells")
            .and_then(|v| v.as_integer())
            .unwrap_or(10) as usize;

        let target_size = min_qc
            .get("target_size")
            .and_then(|v| v.as_real())
            .unwrap_or(1e5) as f32;

        MinCellQuality {
            min_unique_genes,
            min_lib_size,
            min_cells,
            target_size,
        }
    }
}

/// Structure that stores HVG information
///
/// ### Fields
///
/// * `mean` - Mean expression of the gene.
/// * `var` - Detected variance of the gene.
/// * `var_exp` - Expected variance of the gene.
/// * `var_std` - Standardised variance of the gene.
#[derive(Clone, Debug)]
pub struct HvgRes {
    pub mean: Vec<f64>,
    pub var: Vec<f64>,
    pub var_exp: Vec<f64>,
    pub var_std: Vec<f64>,
}

///////////////////////////////
// QC metrics based on genes //
///////////////////////////////

/// Calculates the percentage within the gene set(s)
///
/// Helper function to calculate QC metrics such as mitochondrial proportions,
/// ribosomal proportions, etc.
///
/// ### Params
///
/// * `f_path` - File path to the binarised format that contains the cell-based
///   data
/// * `gene_indices` - Vector of index positions of the genes of interest
/// * `verbose` - Controls verbosity of the function.
///
/// ### Returns
///
/// A vector with the percentages of these genes over the total reads.
pub fn get_gene_set_perc(
    f_path: &str,
    gene_indices: Vec<Vec<u16>>,
    verbose: bool,
) -> Vec<Vec<f32>> {
    let start_reading = Instant::now();

    let reader = ParallelSparseReader::new(f_path).unwrap();

    let cell_chunks = reader.get_all_cells();

    let end_read = start_reading.elapsed();

    if verbose {
        println!("Load in data: {:.2?}", end_read);
    }

    let start_calculations = Instant::now();

    let mut results: Vec<Vec<f32>> = Vec::with_capacity(gene_indices.len());

    for gene_set in gene_indices {
        let hash_gene_set: FxHashSet<&u16> = gene_set.iter().collect();

        let percentage: &Vec<f32> = &cell_chunks
            .par_iter()
            .map(|chunk| {
                let total_sum = chunk
                    .indices
                    .iter()
                    .zip(&chunk.data_raw)
                    .filter(|(col_idx, _)| hash_gene_set.contains(col_idx))
                    .map(|(_, val)| val)
                    .sum::<u16>() as f32;
                let lib_size = chunk.library_size as f32;
                total_sum / lib_size
            })
            .collect();

        results.push(percentage.clone());
    }

    let end_calculations = start_calculations.elapsed();

    if verbose {
        println!(
            "Finished the gene set proportion calculations: {:.2?}",
            end_calculations
        );
    }

    results
}

/// Calculates the percentage within the gene set(s)
///
/// Helper function to calculate QC metrics such as mitochondrial proportions,
/// ribosomal proportions, etc. This function implements streaming and reads in
/// the cells in chunks to avoid memory pressure.
///
/// ### Params
///
/// * `f_path` - File path to the binarised format that contains the cell-based
///   data
/// * `gene_indices` - Vector of index positions of the genes of interest
/// * `verbose` - Controls verbosity of the function.
///
/// ### Returns
///
/// A vector with the percentages of these genes over the total reads.
pub fn get_gene_set_perc_streaming(
    f_path: &str,
    gene_indices: Vec<Vec<u16>>,
    verbose: bool,
) -> Vec<Vec<f32>> {
    let start_total = Instant::now();

    let reader = ParallelSparseReader::new(f_path).unwrap();
    let no_cells = reader.get_header().total_cells;

    const CELL_BATCH_SIZE: usize = 100000;

    let num_cell_batches = no_cells.div_ceil(CELL_BATCH_SIZE);

    let mut results: Vec<Vec<f32>> = vec![Vec::new(); gene_indices.len()];

    let hash_gene_sets: Vec<FxHashSet<&u16>> =
        gene_indices.iter().map(|gs| gs.iter().collect()).collect();

    for cell_batch_idx in 0..num_cell_batches {
        let cell_start = cell_batch_idx * CELL_BATCH_SIZE;
        let cell_end = ((cell_batch_idx + 1) * CELL_BATCH_SIZE).min(no_cells);

        let cell_chunks = reader.read_cells_range(cell_start, cell_end);

        for (gs_idx, hash_gene_set) in hash_gene_sets.iter().enumerate() {
            let percentage: &Vec<f32> = &cell_chunks
                .par_iter()
                .map(|chunk| {
                    let total_sum = chunk
                        .indices
                        .iter()
                        .zip(&chunk.data_raw)
                        .filter(|(col_idx, _)| hash_gene_set.contains(col_idx))
                        .map(|(_, val)| val)
                        .sum::<u16>() as f32;
                    let lib_size = chunk.library_size as f32;
                    total_sum / lib_size
                })
                .collect();

            results[gs_idx].extend(percentage);
        }

        if verbose && cell_batch_idx % 5 == 0 {
            let progress = (cell_batch_idx + 1) as f32 / num_cell_batches as f32 * 100.0;
            println!(
                "  Reading cells and calculating proportions: {:.1}%",
                progress
            );
        }
    }

    let end_total = start_total.elapsed();

    if verbose {
        println!(
            "Finished the gene set proportion calculations: {:.2?}",
            end_total
        );
    }

    results
}

/////////
// HVG //
/////////

/// Calculate the mean and variance of a CSC gene chunk
///
/// ### Params
///
/// * `chunk` - The `CscGeneChunk` representing the gene for which to calculate
///   the mean and variance
/// * `no_cells` - Number of total cells represented in the experiment
///
/// ### Returns
///
/// A tuple of `(mean, var)`
#[inline]
pub fn calculate_mean_var_filtered(
    gene: &CscGeneChunk,
    cell_idx_map: &FxHashMap<u32, u32>,
    no_cells: usize,
) -> (f32, f32) {
    let no_cells = no_cells as f32;
    let mut sum = 0f32;
    let mut nnz = 0usize;

    // Only process cells that are in the filter
    for i in 0..gene.indices.len() {
        if cell_idx_map.contains_key(&gene.indices[i]) {
            sum += gene.data_raw[i] as f32;
            nnz += 1;
        }
    }

    let n_zero = no_cells - nnz as f32;
    let mean = sum / no_cells;

    let mut sum_sq_diff = 0f32;
    for i in 0..gene.indices.len() {
        if cell_idx_map.contains_key(&gene.indices[i]) {
            let val = gene.data_raw[i] as f32;
            let diff = val - mean;
            sum_sq_diff += diff * diff;
        }
    }

    let var = (sum_sq_diff + n_zero * mean * mean) / no_cells;

    (mean, var)
}

/// Helper function to calculate the standardised variance of a gene chunk
///
/// ### Params
///
/// * `chunk` - The `CscGeneChunk` representing the gene for which to calculate
///   the standardised mean and variance
/// * `mean` - Mean value for that gene
/// * `expected_var` - Expected variance based on the Loess function
/// * `clip_max` - Which values to clip
/// * `no_cells` - The number of represented cells
///
/// ### Returns
///
/// The standardised variance
#[inline]
pub fn calculate_std_variance_filtered(
    gene: &CscGeneChunk,
    cell_idx_map: &FxHashMap<u32, u32>,
    mean: f32,
    expected_var: f32,
    clip_max: f32,
    no_cells: usize,
) -> f32 {
    let no_cells_f32 = no_cells as f32;
    let expected_sd = expected_var.sqrt();

    let mut sum_standardised = 0f32;
    let mut sum_sq_standardised = 0f32;
    let mut nnz = 0usize;

    // Process non-zero entries that pass filter
    for i in 0..gene.indices.len() {
        if cell_idx_map.contains_key(&gene.indices[i]) {
            let val_f32 = gene.data_raw[i] as f32;
            let norm = ((val_f32 - mean) / expected_sd)
                .min(clip_max)
                .max(-clip_max);
            sum_standardised += norm;
            sum_sq_standardised += norm * norm;
            nnz += 1;
        }
    }

    // Process zero entries
    let n_zeros = no_cells - nnz;
    if n_zeros > 0 {
        let standardised_zero = ((-mean) / expected_sd).min(clip_max).max(-clip_max);
        sum_standardised += n_zeros as f32 * standardised_zero;
        sum_sq_standardised += n_zeros as f32 * standardised_zero * standardised_zero;
    }

    let standardised_mean = sum_standardised / no_cells_f32;
    (sum_sq_standardised / no_cells_f32) - (standardised_mean * standardised_mean)
}

/// Implementation of the variance stabilised version of the HVG selection
///
/// ### Params
///
/// * `f_path` - Path to the gene-based binary file
/// * `cell_indices` - HashSet with the cell indices to keep.
/// * `loess_span` - Span parameter for the loess function
/// * `clip_max` - Optional clip max parameter
/// * `verbose` - If verbose, returns the timings of the function.
///
/// ### Returns
///
/// The `HvgRes`
pub fn get_hvg_vst(
    f_path: &str,
    cell_indices: &[usize],
    loess_span: f64,
    clip_max: Option<f32>,
    verbose: bool,
) -> HvgRes {
    let start_total = Instant::now();

    // Get data
    let start_read = Instant::now();

    let reader = ParallelSparseReader::new(f_path).unwrap();
    let mut gene_chunks: Vec<CscGeneChunk> = reader.get_all_genes();
    let no_cells = cell_indices.len();

    // build cell mapping ONCE... Before I was doing stupid shit
    let cell_idx_map: FxHashMap<u32, u32> = cell_indices
        .iter()
        .enumerate()
        .map(|(new_idx, &old_idx)| (old_idx as u32, new_idx as u32))
        .collect();

    let end_read = start_read.elapsed();

    if verbose {
        println!("Load in data: {:.2?}", end_read);
    }

    let start_gene_stats = Instant::now();

    let results: Vec<(f32, f32)> = gene_chunks
        .par_iter_mut()
        .map(|chunk| calculate_mean_var_filtered(chunk, &cell_idx_map, no_cells))
        .collect();

    let end_gene_stats = start_gene_stats.elapsed();

    if verbose {
        println!("Calculated gene statistics: {:.2?}", end_gene_stats);
    }

    let start_loess = Instant::now();

    let (means, vars): (Vec<f32>, Vec<f32>) = results.into_iter().unzip();

    let clip_max = clip_max.unwrap_or((no_cells as f32).sqrt());
    let means_log10: Vec<f32> = means.iter().map(|x| x.log10()).collect();
    let vars_log10: Vec<f32> = vars.iter().map(|x| x.log10()).collect();

    let loess = LoessRegression::new(loess_span, 2);
    let loess_res = loess.fit(&means_log10, &vars_log10);

    let end_loess = start_loess.elapsed();

    if verbose {
        println!("Fitted Loess: {:.2?}", end_loess);
    }

    let start_standard = Instant::now();

    let var_standardised: Vec<f32> = gene_chunks
        .par_iter()
        .zip(loess_res.fitted_vals.par_iter())
        .zip(means.par_iter())
        .map(|((chunk_i, var_i), mean_i)| {
            let expected_var = 10_f64.powf(*var_i) as f32;
            calculate_std_variance_filtered(
                chunk_i,
                &cell_idx_map,
                *mean_i,
                expected_var,
                clip_max,
                no_cells,
            )
        })
        .collect();

    let end_standard = start_standard.elapsed();

    if verbose {
        println!("Standardised variance: {:.2?}", end_standard);
    }

    let total = start_total.elapsed();

    if verbose {
        println!("Total run time HVG detection: {:.2?}", total);
    }

    // transform to f64 for R
    HvgRes {
        mean: means.iter().map(|x| *x as f64).collect(),
        var: vars.iter().map(|x| *x as f64).collect(),
        var_exp: loess_res.fitted_vals,
        var_std: var_standardised.iter().map(|x| *x as f64).collect(),
    }
}

/// Implementation of the variance stabilised version of the HVG selection
///
/// This uses a two-pass approach to minimise memory usage:
/// - Pass 1: Calculate mean/variance for loess fitting (genes processed in
///   batches)
/// - Pass 2: Calculate standardized variance using loess results (genes
///   processed in batches)
///
/// ### Params
///
/// * `f_path` - Path to the gene-based binary file
/// * `cell_indices` - Slice with the cell indices to keep.
/// * `loess_span` - Span parameter for the loess function
/// * `clip_max` - Optional clip max parameter
/// * `verbose` - If verbose, returns the timings of the function.
///
/// ### Returns
///
/// The `HvgRes`
pub fn get_hvg_vst_streaming(
    f_path: &str,
    cell_indices: &[usize],
    loess_span: f64,
    clip_max: Option<f32>,
    verbose: bool,
) -> HvgRes {
    let start_total = Instant::now();

    let reader = ParallelSparseReader::new(f_path).unwrap();
    let header = reader.get_header();
    let no_genes = header.total_genes;
    let no_cells = cell_indices.len();

    let cell_idx_map: FxHashMap<u32, u32> = cell_indices
        .iter()
        .enumerate()
        .map(|(new_idx, &old_idx)| (old_idx as u32, new_idx as u32))
        .collect();

    if verbose {
        println!(
            "Pass 1/2: Calculating mean and variance for {} genes...",
            no_genes.separate_with_underscores()
        );
    }

    // Pass 1: Calculate mean and variance in batches
    let start_pass1 = Instant::now();

    const GENE_BATCH_SIZE: usize = 1000;
    let num_batches = no_genes.div_ceil(GENE_BATCH_SIZE);

    let mut means = Vec::with_capacity(no_genes);
    let mut vars = Vec::with_capacity(no_genes);

    for batch_idx in 0..num_batches {
        if verbose && batch_idx % 5 == 0 {
            let progress = (batch_idx + 1) as f32 / num_batches as f32 * 100.0;
            println!("  Progress: {:.1}%", progress);
        }

        let start_gene = batch_idx * GENE_BATCH_SIZE;
        let end_gene = ((batch_idx + 1) * GENE_BATCH_SIZE).min(no_genes);
        let gene_indices: Vec<usize> = (start_gene..end_gene).collect();

        let start_loading = Instant::now();

        let mut genes = reader.read_gene_parallel(&gene_indices);

        let end_loading = start_loading.elapsed();

        if verbose {
            println!("   Loaded batch in: {:.2?}.", end_loading);
        }

        let start_batch = Instant::now();

        let batch_results: Vec<(f32, f32)> = genes
            .par_iter_mut()
            .map(|gene| calculate_mean_var_filtered(gene, &cell_idx_map, no_cells))
            .collect();

        let end_batch = start_batch.elapsed();

        if verbose {
            println!("   Finished calculations in: {:.2?}.", end_batch);
        }

        for (mean, var) in batch_results {
            means.push(mean);
            vars.push(var);
        }
        // genes vec dropped here - memory freed
    }

    let end_pass1 = start_pass1.elapsed();

    if verbose {
        println!("  Calculated gene statistics: {:.2?}", end_pass1);
    }

    // Fit loess
    let start_loess = Instant::now();

    let clip_max = clip_max.unwrap_or((no_cells as f32).sqrt());
    let means_log10: Vec<f32> = means.iter().map(|x| x.log10()).collect();
    let vars_log10: Vec<f32> = vars.iter().map(|x| x.log10()).collect();

    let loess = LoessRegression::new(loess_span, 2);
    let loess_res = loess.fit(&means_log10, &vars_log10);

    let end_loess = start_loess.elapsed();

    if verbose {
        println!("  Fitted Loess: {:.2?}", end_loess);
        println!("Pass 2/2: Calculating standardised variance...");
    }

    // Pass 2: Calculate standardised variance in batches
    let start_pass2 = Instant::now();

    let mut var_standardised = Vec::with_capacity(no_genes);

    for batch_idx in 0..num_batches {
        if verbose && batch_idx % 5 == 0 {
            let progress = (batch_idx + 1) as f32 / num_batches as f32 * 100.0;
            println!("  Progress: {:.1}%", progress);
        }

        let start_gene = batch_idx * GENE_BATCH_SIZE;
        let end_gene = ((batch_idx + 1) * GENE_BATCH_SIZE).min(no_genes);
        let gene_indices: Vec<usize> = (start_gene..end_gene).collect();

        let mut genes = reader.read_gene_parallel(&gene_indices);

        let start_batch = Instant::now();

        let batch_std_vars: Vec<f32> = genes
            .par_iter_mut()
            .enumerate()
            .map(|(local_idx, gene)| {
                let gene_idx = start_gene + local_idx;
                let expected_var = 10_f64.powf(loess_res.fitted_vals[gene_idx]) as f32;
                calculate_std_variance_filtered(
                    gene,
                    &cell_idx_map,
                    means[gene_idx],
                    expected_var,
                    clip_max,
                    no_cells,
                )
            })
            .collect();

        let end_batch = start_batch.elapsed();

        if verbose {
            println!(
                "   Finished calculating standardised variance in: {:.2?}.",
                end_batch
            );
        }

        var_standardised.extend(batch_std_vars);
    }

    let end_pass2 = start_pass2.elapsed();

    if verbose {
        println!(
            "  Calculated standardised variance total: {:.2?}",
            end_pass2
        );
    }

    let total = start_total.elapsed();

    if verbose {
        println!("Total run time HVG detection: {:.2?}", total);
    }

    // transform to f64 for R
    HvgRes {
        mean: means.iter().map(|x| *x as f64).collect(),
        var: vars.iter().map(|x| *x as f64).collect(),
        var_exp: loess_res.fitted_vals,
        var_std: var_standardised.iter().map(|x| *x as f64).collect(),
    }
}

/// To be implemented
pub fn get_hvg_dispersion() -> HvgRes {
    todo!("Dispersion method not yet implemented");

    #[allow(unreachable_code)]
    HvgRes {
        mean: Vec::new(),
        var: Vec::new(),
        var_exp: Vec::new(),
        var_std: Vec::new(),
    }
}

/// To be implemented
pub fn get_hvg_dispersion_streaming() -> HvgRes {
    todo!("Dispersion method with streaming not yet implemented");

    #[allow(unreachable_code)]
    HvgRes {
        mean: Vec::new(),
        var: Vec::new(),
        var_exp: Vec::new(),
        var_std: Vec::new(),
    }
}

/// To be implemented
pub fn get_hvg_mvb() -> HvgRes {
    todo!("MeanVarianceBin method not yet implemented");

    #[allow(unreachable_code)]
    HvgRes {
        mean: Vec::new(),
        var: Vec::new(),
        var_exp: Vec::new(),
        var_std: Vec::new(),
    }
}

/// To be implemented
pub fn get_hvg_mvb_streaming() -> HvgRes {
    todo!("MeanVarianceBin method with streaming not yet implemented");

    #[allow(unreachable_code)]
    HvgRes {
        mean: Vec::new(),
        var: Vec::new(),
        var_exp: Vec::new(),
        var_std: Vec::new(),
    }
}

/////////
// PCA //
/////////

/// Scales the data in a CSC chunk
///
/// ### Params
///
/// * `chunk` - The CscGeneChunk for which to scale the data
/// * `no_cells` - Number of cells represented
///
/// ### Returns
///
/// A densified, scaled vector per gene basis.
#[inline]
fn scale_csc_chunk(chunk: &CscGeneChunk, no_cells: usize) -> Vec<f32> {
    let mut dense_data = vec![0.0f32; no_cells];
    for (idx, &row_idx) in chunk.indices.iter().enumerate() {
        dense_data[row_idx as usize] = chunk.data_norm[idx].to_f32();
    }
    let mean = dense_data.iter().sum::<f32>() / no_cells as f32;
    let variance = dense_data.iter().map(|&x| (x - mean).powi(2)).sum::<f32>() / no_cells as f32;
    let std_dev = variance.sqrt();

    if std_dev < 1e-8 {
        // Zero variance gene - just return centered data (all zeros after centering)
        return vec![0.0f32; no_cells];
    }

    dense_data.iter().map(|&x| (x - mean) / std_dev).collect()
}

/// Calculate the PCs for single cell data
///
/// ### Params
///
/// * `f_path` - Path to the gene-based binary file.
/// * `cell_indices` - HashSet with the cell indices to keep.
/// * `no_pcs` - Number of principal components to calculate
/// * `random_svd` - Shall randomised singular value decompostion be used. This
///   has the advantage of speed-ups, but loses precision.
/// * `return_scaled` - Return the scaled data.
/// * `seed` - Seed for randomised SVD.
///
/// ### Return
///
/// A tuple of the samples projected on the PC space and gene loadings
#[allow(clippy::too_many_arguments)]
pub fn pca_on_sc(
    f_path: &str,
    cell_indices: &[usize],
    gene_indices: &[usize],
    no_pcs: usize,
    random_svd: bool,
    seed: usize,
    return_scaled: bool,
    verbose: bool,
) -> (Mat<f32>, Mat<f32>, Option<Mat<f32>>) {
    let start_total = Instant::now();

    let cell_set: IndexSet<u32> = cell_indices.iter().map(|&x| x as u32).collect();

    let start_reading = Instant::now();

    let reader = ParallelSparseReader::new(f_path).unwrap();
    let mut gene_chunks: Vec<CscGeneChunk> = reader.read_gene_parallel(gene_indices);

    let end_reading = start_reading.elapsed();

    if verbose {
        println!("Loaded in data : {:.2?}", end_reading);
    }

    let start_scaling = Instant::now();

    gene_chunks.par_iter_mut().for_each(|chunk| {
        chunk.filter_selected_cells(&cell_set);
    });

    let scaled_data: Vec<Vec<f32>> = gene_chunks
        .par_iter()
        .map(|chunk| scale_csc_chunk(chunk, cell_indices.len()))
        .collect();

    let num_genes = scaled_data.len();
    let scaled_data = Mat::from_fn(cell_indices.len(), num_genes, |row, col| {
        scaled_data[col][row]
    });

    let end_scaling = start_scaling.elapsed();

    if verbose {
        println!("Finished scaling : {:.2?}", end_scaling);
    }

    let start_svd = Instant::now();

    let (scores, loadings) = if random_svd {
        let res: RandomSvdResults<f32> =
            randomised_svd_f32(scaled_data.as_ref(), no_pcs, seed, Some(100_usize), None);
        // Take first no_pcs components and compute scores as X * V
        let loadings = res.v.submatrix(0, 0, num_genes, no_pcs).to_owned();
        let scores = &scaled_data * &loadings;
        (scores, loadings)
    } else {
        let res = scaled_data.thin_svd().unwrap();
        // Take only the first no_pcs components
        let loadings = res.V().submatrix(0, 0, num_genes, no_pcs).to_owned();
        let scores = &scaled_data * &loadings;
        (scores, loadings)
    };

    let end_svd = start_svd.elapsed();

    if verbose {
        println!("Finished PCA calculations : {:.2?}", end_svd);
    }

    let end_total = start_total.elapsed();

    if verbose {
        println!("Total run time PCA detection: {:.2?}", end_total);
    }

    let scaled = if return_scaled {
        Some(scaled_data)
    } else {
        None
    };

    (scores, loadings, scaled)
}
