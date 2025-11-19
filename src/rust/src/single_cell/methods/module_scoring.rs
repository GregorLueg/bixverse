use indexmap::IndexSet;
use rand::prelude::IndexedRandom;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::time::Instant;

use crate::core::data::sparse_io::*;

/////////////
// Helpers //
/////////////

/// Structure holding gene binning information
///
/// ### Fields
///
/// * `gene_to_bin` - HashMap that contains gene to bin mapping
/// * `bins` - Bin to gene lookup
struct GeneBins {
    gene_to_bin: FxHashMap<usize, usize>,
    bins: Vec<Vec<usize>>,
}

/// Helper function to get the average gene expression across the cells
///
/// ### Params
///
/// * `f_path_gene` - Path to the gene-based binary file.
/// * `cell_set` - IndexSet that stores which cells to include in the analysis.
/// * `streaming` - Boolean. If set to TRUE, the chunks will be loaded in groups
///   of 500 gene.
///
/// ### Returns
///
/// A vector of `(gene_index, avg expression)`
fn get_average_expression(
    f_path_gene: &str,
    cell_set: &IndexSet<u32>,
    streaming: bool,
) -> Result<Vec<(usize, f32)>, String> {
    let reader = ParallelSparseReader::new(f_path_gene)
        .map_err(|e| format!("Failed to open file: {}", e))?;
    let total_genes = reader.get_header().total_genes;

    if streaming {
        // 500 genes in one go
        const CHUNK_SIZE: usize = 500;
        let gene_indices: Vec<usize> = (0..total_genes).collect();

        let results: Vec<(usize, f32)> = gene_indices
            .chunks(CHUNK_SIZE)
            .flat_map(|chunk| {
                let gene_chunks = reader.read_gene_parallel(chunk);
                gene_chunks
                    .par_iter()
                    .map(|gene| gene.calculate_avg_exp(cell_set))
                    .collect::<Vec<_>>()
            })
            .collect();

        Ok(results)
    } else {
        let gene_chunks = reader.get_all_genes();
        let results: Vec<(usize, f32)> = gene_chunks
            .par_iter()
            .map(|gene| gene.calculate_avg_exp(cell_set))
            .collect();

        Ok(results)
    }
}

/// Create expression bins for genes
///
/// Bins genes into equal-sized groups based on average expression.
/// Follows Seurat's approach using quantile-based binning.
///
/// ### Params
///
/// * `gene_means` - Slice of tuples (gene_index, average_expression)
/// * `nbin` - Number of bins to create
///
/// ### Returns
///
/// `GeneBins` structure with gene->bin mapping and bin->genes lookup
fn create_expression_bins(gene_means: &[(usize, f32)], nbin: usize, seed: &usize) -> GeneBins {
    let mut sorted_genes = gene_means.to_vec();

    // add tiny random noise to break ties
    // seurat does this
    let mut rng = StdRng::seed_from_u64(*seed as u64);
    for (_, exp) in sorted_genes.iter_mut() {
        *exp += rng.random::<f32>() / 1e30;
    }

    // Sort by expression
    sorted_genes.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

    // Create equal-sized bins
    let total_genes = sorted_genes.len();
    let genes_per_bin = (total_genes as f32 / nbin as f32).ceil() as usize;

    let mut gene_to_bin = FxHashMap::default();
    let mut bins: Vec<Vec<usize>> = vec![Vec::new(); nbin];

    for (i, (gene_idx, _)) in sorted_genes.iter().enumerate() {
        let bin_id = (i / genes_per_bin).min(nbin - 1);
        gene_to_bin.insert(*gene_idx, bin_id);
        bins[bin_id].push(*gene_idx);
    }

    GeneBins { gene_to_bin, bins }
}

/// Sample control genes for a gene set
///
/// For each gene in the set, samples `ctrl` genes from the same expression bin.
///
/// ### Params
///
/// * `gene_set` - Slice of gene indices in the set
/// * `gene_bins` - Gene binning structure
/// * `ctrl` - Number of control genes per feature
/// * `rng` - Random number generator
///
/// ### Returns
///
/// Vec of unique control gene indices
fn sample_control_genes(
    gene_set: &[usize],
    gene_bins: &GeneBins,
    ctrl: usize,
    rng: &mut impl Rng,
) -> Vec<usize> {
    let mut controls = Vec::new();

    for &gene_idx in gene_set {
        if let Some(&bin_id) = gene_bins.gene_to_bin.get(&gene_idx) {
            let bin_genes = &gene_bins.bins[bin_id];

            // Sample ctrl genes from this bin (with replacement as per Seurat)
            let sampled: Vec<usize> = bin_genes
                .choose_multiple(rng, ctrl.min(bin_genes.len()))
                .copied()
                .collect();
            controls.extend(sampled);
        }
    }

    // Return unique controls
    controls.sort_unstable();
    controls.dedup();
    controls
}

/// Calculate module scores for a single cell
///
/// ### Params
///
/// * `cell` - Reference to a CsrCellChunk
/// * `gene_set` - Indices of genes in the module
/// * `control_set` - Indices of control genes
///
/// ### Returns
///
/// Module score defined as `(mean(genes_of_interest) - mean(controls))`
fn calculate_cell_module_score(
    cell: &CsrCellChunk,
    gene_set: &[usize],
    control_set: &[usize],
) -> f32 {
    let find_expr = |gene_idx: usize| -> f32 {
        cell.indices
            .binary_search(&(gene_idx as u16))
            .ok()
            .map(|pos| cell.data_norm[pos].to_f32())
            .unwrap_or(0.0)
    };

    let gene_sum: f32 = gene_set.iter().map(|&idx| find_expr(idx)).sum();
    let gene_mean = if gene_set.is_empty() {
        0.0
    } else {
        gene_sum / gene_set.len() as f32
    };

    let ctrl_sum: f32 = control_set.iter().map(|&idx| find_expr(idx)).sum();
    let ctrl_mean = if control_set.is_empty() {
        0.0
    } else {
        ctrl_sum / control_set.len() as f32
    };

    gene_mean - ctrl_mean
}

/// Calculate the module scores
///
/// ### Params
///
/// * `f_path_cell` - Path to the cell-based binary file
/// * `gene_sets` - Slice of indices of the gene sets
/// * `cells_to_keep` - Slice of indices of the cells to keep
/// * `gene_bins` - The pre-calculated GeneBins.
/// * `ctrl` - Number of control genes to use
/// * `seed` - Seed for reproducibility
///
/// ### Returns
///
/// Vec of vec with outer vector representing the gene sets and the inner ones
/// the cells.
fn calculate_module_scores(
    f_path_cell: &str,
    gene_sets: &[Vec<usize>],
    cells_to_keep: &[usize],
    gene_bins: &GeneBins,
    ctrl: usize,
    seed: &usize,
) -> Result<Vec<Vec<f32>>, String> {
    let mut rng = StdRng::seed_from_u64(*seed as u64);

    let control_sets: Vec<Vec<usize>> = gene_sets
        .iter()
        .map(|gene_set| sample_control_genes(gene_set, gene_bins, ctrl, &mut rng))
        .collect();

    let reader = ParallelSparseReader::new(f_path_cell)
        .map_err(|e| format!("Failed to open file: {}", e))?;
    let cell_chunks = reader.read_cells_parallel(cells_to_keep);

    let all_scores: Vec<Vec<f32>> = cell_chunks
        .par_iter()
        .map(|cell| {
            gene_sets
                .iter()
                .zip(control_sets.iter())
                .map(|(gene_set, control_set)| {
                    calculate_cell_module_score(cell, gene_set, control_set)
                })
                .collect()
        })
        .collect();

    // Transpose: cells x modules -> modules x cells
    let mut results: Vec<Vec<f32>> = vec![Vec::with_capacity(cell_chunks.len()); gene_sets.len()];
    for cell_scores in all_scores {
        for (module_idx, score) in cell_scores.into_iter().enumerate() {
            results[module_idx].push(score);
        }
    }

    Ok(results)
}

/// Calculate the module scores (in a streaming fashion)
///
/// ### Params
///
/// * `f_path_cell` - Path to the cell-based binary file
/// * `gene_sets` - Slice of indices of the gene sets
/// * `cells_to_keep` - Slice of indices of the cells to keep
/// * `gene_bins` - The pre-calculated GeneBins.
/// * `ctrl` - Number of control genes to use
/// * `seed` - Seed for reproducibility
///
/// ### Returns
///
/// Vec of vec with outer vector representing the gene sets and the inner ones
/// the cells.
fn calculate_module_scores_streaming(
    f_path_cell: &str,
    gene_sets: &[Vec<usize>],
    cells_to_keep: &[usize],
    gene_bins: &GeneBins,
    ctrl: usize,
    seed: &usize,
    verbose: bool,
) -> Result<Vec<Vec<f32>>, String> {
    const CHUNK_SIZE: usize = 50000;

    let mut rng = StdRng::seed_from_u64(*seed as u64);

    let control_sets: Vec<Vec<usize>> = gene_sets
        .iter()
        .map(|gene_set| sample_control_genes(gene_set, gene_bins, ctrl, &mut rng))
        .collect();

    let reader = ParallelSparseReader::new(f_path_cell)
        .map_err(|e| format!("Failed to open file: {}", e))?;

    let total_chunks = cells_to_keep.len().div_ceil(CHUNK_SIZE);
    let mut results: Vec<Vec<f32>> = vec![Vec::with_capacity(cells_to_keep.len()); gene_sets.len()];

    for (chunk_idx, cell_indices_chunk) in cells_to_keep.chunks(CHUNK_SIZE).enumerate() {
        let start = Instant::now();

        let cell_chunks = reader.read_cells_parallel(cell_indices_chunk);

        // Calculate scores in parallel (cells x modules)
        let chunk_scores: Vec<Vec<f32>> = cell_chunks
            .par_iter()
            .map(|cell| {
                gene_sets
                    .iter()
                    .zip(control_sets.iter())
                    .map(|(gene_set, control_set)| {
                        calculate_cell_module_score(cell, gene_set, control_set)
                    })
                    .collect()
            })
            .collect();

        // Transpose and append: cells x modules -> modules x cells
        for cell_scores in chunk_scores {
            for (module_idx, score) in cell_scores.into_iter().enumerate() {
                results[module_idx].push(score);
            }
        }

        if verbose {
            let elapsed = start.elapsed();
            let pct = ((chunk_idx + 1) as f32 / total_chunks as f32) * 100.0;
            println!(
                "Chunk {} of {} (took {:.2?}, {:.1}% complete)",
                chunk_idx + 1,
                total_chunks,
                elapsed,
                pct
            );
        }
    }

    Ok(results)
}

/// Calculate the module scores
///
/// ### Params
///
/// * `f_path_gene` - Path to the gene-based binary file.
/// * `f_path_cell` - Path to the cell-based binary file.
/// * `gene_sets` - Slice of indices of the gene sets.
/// * `cells_to_use` - Slice of indices of the cells to use.
/// * `nbin` - Number of bins to use
/// * `ctrl` - Number of control genes to use.
/// * `streaming` - Shall streaming be used. Useful for larger data sets.
/// * `seed` - Seed for reproducibility.
/// * `verbose` - Controls verbosity of the function.
///
/// ### Returns
///
/// Vec of vec with outer vector representing the gene sets and the inner ones
/// the cells.
#[allow(clippy::too_many_arguments)]
pub fn calculate_module_scores_main(
    f_path_gene: &str,
    f_path_cell: &str,
    gene_sets: &[Vec<usize>],
    cells_to_use: &[usize],
    nbin: usize,
    ctrl: usize,
    streaming: bool,
    seed: usize,
    verbose: bool,
) -> Result<Vec<Vec<f32>>, String> {
    let cell_set: IndexSet<u32> = cells_to_use.iter().map(|&x| x as u32).collect();

    let start_total = Instant::now();

    if verbose {
        println!("Calculating the average expression across the cells.")
    }

    let start_avg_exp = Instant::now();

    let avg_exp = get_average_expression(f_path_gene, &cell_set, streaming)?;

    let end_evg_exp = start_avg_exp.elapsed();

    if verbose {
        println!(
            "Finished the calculation of the avg gene expression in {:.2?}",
            end_evg_exp
        );
        println!("Calculating the module scores now.")
    }

    let start_modules = Instant::now();

    let gene_bins = create_expression_bins(&avg_exp, nbin, &seed);

    let module_scores = if streaming {
        calculate_module_scores_streaming(
            f_path_cell,
            gene_sets,
            cells_to_use,
            &gene_bins,
            ctrl,
            &seed,
            verbose,
        )
    } else {
        calculate_module_scores(
            f_path_cell,
            gene_sets,
            cells_to_use,
            &gene_bins,
            ctrl,
            &seed,
        )
    }?;

    let end_modules = start_modules.elapsed();
    let end_total = start_total.elapsed();

    if verbose {
        println!(
            "Finished the calculation of the modules in {:.2?}",
            end_modules
        );
        println!("Total runtime: {:.2?}", end_total)
    }

    Ok(module_scores)
}
