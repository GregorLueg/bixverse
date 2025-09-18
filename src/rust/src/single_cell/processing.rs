use extendr_api::*;
use rayon::prelude::*;
use rustc_hash::FxHashSet;

use crate::core::base::loess::*;
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
///
/// ### Returns
///
/// A vector with the percentages of these genes over the total reads.
pub fn get_gene_set_perc(f_path: &str, gene_indices: Vec<Vec<u16>>) -> Vec<Vec<(f32, f32, f32)>> {
    let reader = ParallelSparseReader::new(f_path).unwrap();

    let cell_chunks = reader.get_all_cells();

    let mut results: Vec<Vec<(f32, f32, f32)>> = Vec::with_capacity(gene_indices.len());

    for gene_set in gene_indices {
        let hash_gene_set: FxHashSet<&u16> = gene_set.iter().collect();

        let percentage: &Vec<(f32, f32, f32)> = &cell_chunks
            .par_iter()
            .map(|chunk| {
                let total_sum = chunk
                    .col_indices
                    .iter()
                    .zip(&chunk.data_raw)
                    .filter(|(col_idx, _)| hash_gene_set.contains(col_idx))
                    .map(|(_, val)| val)
                    .sum::<u16>() as f32;
                let lib_size = chunk.library_size as f32;
                (total_sum / lib_size, total_sum, lib_size)
            })
            .collect();

        results.push(percentage.clone());
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
fn calculate_mean_var_csc_chunk(chunk: &CscGeneChunk, no_cells: usize) -> (f32, f32) {
    let no_cells = no_cells as f32;
    let nnz = chunk.nnz as f32;
    let sum: f32 = chunk.data_raw.iter().map(|x| *x as f32).sum();
    let n_zero = no_cells - nnz;

    let mean = sum / no_cells;
    let sum_sq_diff = chunk
        .data_raw
        .iter()
        .map(|&val| {
            let val: f32 = val.into();
            let diff = val - mean;
            diff * diff
        })
        .sum::<f32>();

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
#[inline]
fn calculate_std_variance(
    chunk: &CscGeneChunk,
    mean: &f32,
    expected_var: &f32,
    clip_max: &f32,
    no_cells: usize,
) -> f32 {
    let no_cells_f32 = no_cells as f32;
    let expected_sd = expected_var.sqrt();

    let mut sum_standardised = 0_f32;
    let mut sum_sq_standardised = 0_f32;

    // process the non-zero entries
    for val in &chunk.data_raw {
        let val_f32 = *val as f32;
        let norm = ((val_f32 - mean) / expected_sd)
            .min(*clip_max)
            .max(-clip_max);
        sum_standardised += norm;
        sum_sq_standardised += norm * norm;
    }

    // process the zero entries
    let n_zeros = no_cells - chunk.nnz;
    if n_zeros > 0 {
        let standardised_zero = ((-mean) / expected_sd).min(*clip_max).max(-clip_max);
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
///
/// ### Returns
///
/// The `HvgRes`
pub fn get_hvg_vst(
    f_path: &str,
    cell_indices: &FxHashSet<u32>,
    loess_span: f64,
    clip_max: Option<f32>,
) -> HvgRes {
    // Get data
    let reader = ParallelSparseReader::new(f_path).unwrap();
    let mut gene_chunks: Vec<CscGeneChunk> = reader.get_all_genes();
    let no_cells = cell_indices.len();

    let results: Vec<(f32, f32)> = gene_chunks
        .par_iter_mut()
        .map(|chunk| {
            chunk.filter_selected_cells(cell_indices);
            calculate_mean_var_csc_chunk(chunk, no_cells)
        })
        .collect();

    let (means, vars): (Vec<f32>, Vec<f32>) = results.into_iter().unzip();

    let clip_max = clip_max.unwrap_or((no_cells as f32).sqrt());
    let means_log10: Vec<f32> = means.iter().map(|x| x.log10()).collect();
    let vars_log10: Vec<f32> = vars.iter().map(|x| x.log10()).collect();

    let loess = LoessRegression::new(loess_span, 1);
    let loess_res = loess.fit(&means_log10, &vars_log10);

    let var_standardised: Vec<f32> = gene_chunks
        .par_iter()
        .zip(loess_res.fitted_vals.par_iter())
        .zip(means.par_iter())
        .map(|((chunk_i, var_i), mean_i)| {
            let var_32 = *var_i as f32;
            calculate_std_variance(chunk_i, mean_i, &var_32, &clip_max, no_cells)
        })
        .collect();

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

    HvgRes {
        mean: Vec::new(),
        var: Vec::new(),
        var_exp: Vec::new(),
        var_std: Vec::new(),
    }
}
