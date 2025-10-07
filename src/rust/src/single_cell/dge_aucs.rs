use faer::Mat;
use half::f16;
use rand::prelude::IndexedRandom;
use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use std::time::Instant;

use crate::core::base::stats::{calc_fdr, z_scores_to_pval};
use crate::core::data::sparse_io::*;
use crate::core::data::sparse_structures::CompressedSparseData;
use crate::single_cell::fast_ranking::rank_csr_chunk_vec;

///////////
// Enums //
///////////

/// Enum describing the type of AUC to calculate
#[derive(Clone, Debug)]
enum AucType {
    /// AUC based on the Mann-Whitney statistic
    MannWhitney,
    /// AUC based on the traditional AUC/AUROC calculation
    ClassicalAuc,
}

////////////////
// Structures //
////////////////

/// Structure to store the Mann Whitney U-based DGE results in
///
/// ### Fields
///
/// * `lfc` - The log folc changes
/// * `prop1` - Proportions of cells in group 1 expressing the gene.
/// * `prop2` - Proportions of cells in group 2 expressing the gene.
/// * `z_scores` - The Z-scores based on the Mann-Whitney U test.
/// * `p_vals` - The p-values from the Mann-Whitney U test.
/// * `fdr` - FDR values given the p-values.
/// * `genes_to_keep` - Boolean indicating if the gene was included in the
///   analysis, i.e., passed the proportion thresholds.
#[derive(Clone, Debug)]
pub struct DgeMannWhitneyRes {
    pub lfc: Vec<f32>,
    pub prop1: Vec<f32>,
    pub prop2: Vec<f32>,
    pub z_scores: Vec<f64>,
    pub p_vals: Vec<f64>,
    pub fdr: Vec<f64>,
    pub genes_to_keep: Vec<bool>,
}

/////////////
// Helpers //
/////////////

/// Calculate the average expression and proportion for the genes
///
/// Takes in slice of CsrCellChunks and calculates the average expression
/// across the genes and proportions in which the gene is expressed.
///
/// ### Params
///
/// * `cells` - A vector of `CsrCellChunk`.
/// * `num_genes` - Number of represented genes in the data.
///
/// ### Returns
///
/// A tuple of `(average expression, proportion of cells expressing gene)`
fn calculate_avg_exp_prop(cells: &[CsrCellChunk], num_genes: usize) -> (Vec<f32>, Vec<f32>) {
    let mut sum_exp = vec![0.0f32; num_genes];
    let mut count_exp = vec![0usize; num_genes];

    for cell in cells {
        for (&gene_idx, &norm_val) in cell.indices.iter().zip(cell.data_norm.iter()) {
            sum_exp[gene_idx as usize] += norm_val.to_f32();
            count_exp[gene_idx as usize] += 1;
        }
    }

    let total_cells = cells.len() as f32;
    let avg_exp: Vec<f32> = sum_exp.iter().map(|&sum| sum / total_cells).collect();
    let prop_exp: Vec<f32> = count_exp
        .iter()
        .map(|&count| count as f32 / total_cells)
        .collect();

    (avg_exp, prop_exp)
}

/// Calculates the Mann Whitney statistic
///
/// ### Params
///
/// * `ranks1` - Ranked data for the gene of group 1
/// * `ranks2` - Ranked data for the gene of group 2
///
/// ### Returns
///
/// Returns the Z-score based on Mann-Whitney test
fn mann_whitney_u_test(ranks1: &[f32], ranks2: &[f32]) -> f64 {
    let n1 = ranks1.len() as f64;
    let n2 = ranks2.len() as f64;

    if n1 == 0.0 || n2 == 0.0 {
        return 0.0;
    }

    let r1: f64 = ranks1.iter().map(|&r| r as f64).sum();
    let u1 = n1 * n2 + n1 * (n1 + 1.0) / 2.0 - r1;

    let mean = n1 * n2 / 2.0;
    let variance = n1 * n2 * (n1 + n2 + 1.0) / 12.0;

    (mean - u1) / variance.sqrt()
}

/// Parse the desired AUC type
///
/// ### Params
///
/// * `s` String specifying the desired AUC type.
///
/// ### Return
///
/// The Option of the `AucType`
fn parse_auc_type(s: &str) -> Option<AucType> {
    match s.to_lowercase().as_str() {
        "auroc" => Some(AucType::ClassicalAuc),
        "wilcox" => Some(AucType::MannWhitney),
        _ => None,
    }
}

/// Calculate AUC based on ranks and gene set indices (Mann-Whitney version)
///
/// This uses the Mann Whitney statistic under the hood and calculates how
/// active the gene set is over a random gene set. Question asked:
/// "Do genes in my set rank higher than genes not in my set?"
///
/// ### Params
///
/// * `ranks` - The within cell ranked data.
/// * `gene_set` - Indices of the members of this gene set.
///
/// ### Returns
///
/// AUC for this gene set based on the Mann Whitney statistc.
fn calculate_auc_per_cell_mw(ranks: &[f32], gene_set: &[usize]) -> f32 {
    let n_genes = ranks.len();
    let n_in_set = gene_set.len();
    let n_not_in_set = n_genes - n_in_set;

    // Sum of ranks for genes in the set
    let rank_sum: f32 = gene_set.iter().map(|&idx| ranks[idx]).sum();

    // U statistic
    let u = rank_sum - (n_in_set * (n_in_set + 1)) as f32 / 2.0;

    // Convert to AUC
    u / (n_in_set * n_not_in_set) as f32
}

/// Calculate AUC based on ranks and gene set indices (classical AUC)
///
/// This function uses a more classical appraoch of the AUC calculations and
/// is more top-heavy sensitive and can be for example used for marker gene
/// detection. Question asked: "How enriched is my gene set at the top of the
/// ranking?"
///
/// ### Params
///
/// * `ranks` - The within cell ranked data.
/// * `gene_set` - Indices of the members of this gene set.
///
/// ### Returns
///
/// AUC for this gene set based on the Mann Whitney statistc.
fn calculate_auc_for_cell_auroc(ranks: &[f32], gene_set: &[usize]) -> f32 {
    let n_genes = ranks.len();
    let mut gene_set_ranks: Vec<f32> = gene_set.iter().map(|&idx| ranks[idx]).collect();

    gene_set_ranks.sort_by(|a, b| a.partial_cmp(b).unwrap());

    // Calculate AUC: sum of (rank_position - expected_position)
    let n_genes_in_set = gene_set_ranks.len() as f32;
    let auc: f32 = gene_set_ranks
        .iter()
        .enumerate()
        .map(|(i, &rank)| {
            let expected_rank = (i + 1) as f32 * n_genes as f32 / n_genes_in_set;
            rank - expected_rank
        })
        .sum();

    // Normalize
    let max_auc = n_genes_in_set * (n_genes as f32 - n_genes_in_set / 2.0);
    auc / max_auc
}

////////////////////
// Main functions //
////////////////////

/// Get differential expression based on Mann-Whitney
///
/// ### Params
///
/// * `f_path` - File path to the cell-based binary file.
/// * `grp_1_indices` - The cell indices of group 1.
/// * `grp_2_indices` - The cell indices of group 2.
/// * `min_proportion` - The minimum proportion that a gene needs to be
///   expressed in at least one of the two groups.
/// * `alternative` - The test alternative. One of `"twosided"`, `"greater"`,
///   or `"less"`
/// * `verbose` - Boolean. Controls the verbosity of the function.
///
/// ### Returns
///
/// The `DgeMannWhitneyRes` structure with results
pub fn calculate_dge_grps_mann_whitney(
    f_path: &str,
    grp_1_indices: &[usize],
    grp_2_indices: &[usize],
    min_proportion: f32,
    alternative: &str,
    verbose: bool,
) -> Result<DgeMannWhitneyRes, String> {
    let start_read = Instant::now();

    let reader = ParallelSparseReader::new(f_path).unwrap();
    let no_genes = reader.get_header().total_genes;

    let mut cell_chunks_1: Vec<CsrCellChunk> = reader.read_cells_parallel(grp_1_indices);
    let mut cell_chunks_2: Vec<CsrCellChunk> = reader.read_cells_parallel(grp_2_indices);

    let end_read = start_read.elapsed();

    if verbose {
        println!("Loaded in data: {:.2?}", end_read);
    }

    let (avg_exp_1, prop_1) = calculate_avg_exp_prop(&cell_chunks_1, no_genes);
    let (avg_exp_2, prop_2) = calculate_avg_exp_prop(&cell_chunks_2, no_genes);

    let genes_to_keep: Vec<bool> = prop_1
        .iter()
        .zip(prop_2.iter())
        .map(|(&p1, &p2)| p1 >= min_proportion || p2 >= min_proportion)
        .collect();

    let no_genes_new = genes_to_keep.iter().filter(|&&x| x).count();

    cell_chunks_1
        .par_iter_mut()
        .for_each(|cell| cell.filter_genes(&genes_to_keep));
    cell_chunks_2
        .par_iter_mut()
        .for_each(|cell| cell.filter_genes(&genes_to_keep));

    let genes_kept: Vec<usize> = genes_to_keep
        .iter()
        .enumerate()
        .filter_map(|(i, &keep)| if keep { Some(i) } else { None })
        .collect();

    let n1 = cell_chunks_1.len();

    let mut combined_chunks = cell_chunks_1;
    combined_chunks.extend(cell_chunks_2);

    let start_ranking = Instant::now();

    let all_ranks = rank_csr_chunk_vec(combined_chunks, no_genes_new, false);

    let end_ranking = start_ranking.elapsed();

    if verbose {
        println!("Finished the ranking across cells: {:.2?}", end_ranking);
    }

    let start_calculations = Instant::now();

    let res: Vec<(f32, f32, f32, f64)> = genes_kept
        .par_iter()
        .enumerate()
        .map(|(new_idx, &original_idx)| {
            let log_fc = avg_exp_1[original_idx] - avg_exp_2[original_idx];
            let prop1 = prop_1[original_idx];
            let prop2 = prop_2[original_idx];
            let gene_ranks = &all_ranks[new_idx];

            let z = mann_whitney_u_test(&gene_ranks[..n1], &gene_ranks[n1..]);

            (log_fc, prop1, prop2, z)
        })
        .collect();

    let mut log_fc = Vec::with_capacity(res.len());
    let mut prop1 = Vec::with_capacity(res.len());
    let mut prop2 = Vec::with_capacity(res.len());
    let mut z_scores = Vec::with_capacity(res.len());

    for (log_fc_i, prop1_i, prop2_i, z_i) in res {
        log_fc.push(log_fc_i);
        prop1.push(prop1_i);
        prop2.push(prop2_i);
        z_scores.push(z_i);
    }

    let p_vals = z_scores_to_pval(&z_scores, alternative)?;
    let fdr = calc_fdr(&p_vals);

    let end_calculations = start_calculations.elapsed();

    if verbose {
        println!("Finished DGE calculations: {:.2?}", end_calculations);
    }

    Ok(DgeMannWhitneyRes {
        lfc: log_fc,
        prop1,
        prop2,
        z_scores,
        p_vals,
        fdr,
        genes_to_keep,
    })
}

/// Calculate AUCell
///
/// Function to calculate AUC values for gene sets on a per gene set basis.
/// Has the option to calculate the AUCs in form of AUROCs (top heavy, useful
/// for marker gene expression) or derived from the Mann Whitney statistic
/// (is the gene set more active compared to the rest of the GEX; useful for
/// pathway activity measures).
///
/// ### Params
///
/// * `f_path` -  File path to the cell-based binary file.
/// * `gene_sets` - Slice of Vecs indicating the indices of the gene sets
/// * `auc_type` - String. One of `"auroc"` or `"wilcox`
///
/// ### Returns
///
/// AUCell-type values in form gene set x cells.
pub fn calculate_aucell(
    f_path: &str,
    gene_sets: &[Vec<usize>],
    auc_type: &str,
    verbose: bool,
) -> Result<Vec<Vec<f32>>, String> {
    let auc_type = parse_auc_type(auc_type)
        .ok_or_else(|| format!("Invalid AUC method: {}", auc_type))
        .unwrap();

    let start_read = Instant::now();
    let reader = ParallelSparseReader::new(f_path).unwrap();
    let no_genes = reader.get_header().total_genes;
    let cell_chunks: Vec<CsrCellChunk> = reader.get_all_cells();
    let total_cells = cell_chunks.len();
    let end_read = start_read.elapsed();

    if verbose {
        println!("Loaded in data: {:.2?}", end_read);
    }

    let start_ranking = Instant::now();
    let ranks = rank_csr_chunk_vec(cell_chunks, no_genes, true);
    let end_ranking = start_ranking.elapsed();

    if verbose {
        println!("Ranked gene expression within cells {:.2?}", end_ranking);
    }

    let start_auc = Instant::now();
    let mut all_results: Vec<Vec<f32>> = vec![Vec::with_capacity(total_cells); gene_sets.len()];

    for cell_ranks in ranks {
        let aucs: Vec<f32> = gene_sets
            .par_iter()
            .map(|gene_set| match auc_type {
                AucType::ClassicalAuc => calculate_auc_for_cell_auroc(&cell_ranks, gene_set),
                AucType::MannWhitney => calculate_auc_per_cell_mw(&cell_ranks, gene_set),
            })
            .collect();

        for (gene_set_idx, auc) in aucs.into_iter().enumerate() {
            all_results[gene_set_idx].push(auc);
        }
    }
    let end_auc = start_auc.elapsed();

    if verbose {
        println!("Calulated AUCs {:.2?}", end_auc);
    }

    Ok(all_results)
}

/// Calculate AUCell (streaming)
///
/// Function to calculate AUC values for gene sets on a per gene set basis.
/// Has the option to calculate the AUCs in form of AUROCs (top heavy, useful
/// for marker gene expression) or derived from the Mann Whitney statistic
/// (is the gene set more active compared to the rest of the GEX; useful for
/// pathway activity measures). This version streams the data in chunks of
/// 50,000 cells to reduce memory pressure.
///
/// ### Params
///
/// * `f_path` -  File path to the cell-based binary file.
/// * `gene_sets` - Slice of Vecs indicating the indices of the gene sets
/// * `auc_type` - String. One of `"auroc"` or `"wilcox`
///
/// ### Returns
///
/// AUCell-type values in form gene set x cells.
pub fn calculate_aucell_streaming(
    f_path: &str,
    gene_sets: &[Vec<usize>],
    auc_type: &str,
    verbose: bool,
) -> Result<Vec<Vec<f32>>, String> {
    const CHUNK_SIZE: usize = 50000;

    let auc_type = parse_auc_type(auc_type)
        .ok_or_else(|| format!("Invalid AUC method: {}", auc_type))
        .unwrap();

    let reader = ParallelSparseReader::new(f_path).unwrap();
    let no_genes = reader.get_header().total_genes;
    let total_cells = reader.get_header().total_cells;
    let total_chunks = total_cells.div_ceil(CHUNK_SIZE);

    let mut all_results: Vec<Vec<f32>> = vec![Vec::with_capacity(total_cells); gene_sets.len()];

    for (chunk_idx, chunk_start) in (0..total_cells).step_by(CHUNK_SIZE).enumerate() {
        let start_chunk = Instant::now();

        let chunk_end = (chunk_start + CHUNK_SIZE).min(total_cells);
        let cell_indices: Vec<usize> = (chunk_start..chunk_end).collect();

        let cell_chunks = reader.read_cells_parallel(&cell_indices);
        let ranks = rank_csr_chunk_vec(cell_chunks, no_genes, true);

        for cell_ranks in ranks {
            for (gene_set_idx, gene_set) in gene_sets.iter().enumerate() {
                let auc = match auc_type {
                    AucType::ClassicalAuc => calculate_auc_for_cell_auroc(&cell_ranks, gene_set),
                    AucType::MannWhitney => calculate_auc_per_cell_mw(&cell_ranks, gene_set),
                };
                all_results[gene_set_idx].push(auc);
            }
        }

        if verbose {
            let elapsed = start_chunk.elapsed();
            let pct_complete = ((chunk_idx + 1) as f32 / total_chunks as f32) * 100.0;
            println!(
                "Processing chunk {} out of {} (took {:.2?}, completed {:.1}%)",
                chunk_idx + 1,
                total_chunks,
                elapsed,
                pct_complete
            );
        }
    }

    Ok(all_results)
}

///////////////////////
// Cell aggregations //
///////////////////////

/////////////
// Helpers //
/////////////

/// Select meta cells
///
/// ### Params
///
/// * `nn_map` - Nearest neighbours with self
/// * `max_shared` - Maximum number of shared neighbours to be considered
/// * `target_no` - Target number of meta cells
/// * `max_iter` - Maximum iterations for the algorithm
/// * `seed` - seed for reproducibility purposes
///
/// ### Returns
///
/// Borrowed slice of the selected meta cell indices
pub fn identify_meta_cells(
    nn_map: &[Vec<usize>],
    max_shared: usize,
    target_no: usize,
    max_iter: usize,
    seed: usize,
) -> Vec<&[usize]> {
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed as u64);
    let k = nn_map[0].len();
    let k2 = k * 2;

    let mut good_choices: Vec<usize> = (0..nn_map.len()).collect();
    let mut chosen: Vec<usize> = Vec::new();

    if let Some(&first) = good_choices.choose(&mut rng) {
        chosen.push(first);
        good_choices.retain(|&x| x != first);
    }

    let mut it = 0;

    // cache the HashSets to avoid regeneration during loops
    let mut set_cache: FxHashMap<usize, FxHashSet<usize>> = FxHashMap::default();

    // bootstrap meta cells
    while !good_choices.is_empty() && chosen.len() < target_no && it < max_iter {
        it += 1;

        // sample remaining cells
        let choice_idx = rng.random_range(0..good_choices.len());
        let candidate = good_choices[choice_idx];
        good_choices.remove(choice_idx);

        // check overlap with existing meta cells
        set_cache
            .entry(candidate)
            .or_insert_with(|| nn_map[candidate].iter().copied().collect());

        let mut max_overlap = 0;
        for &existing in &chosen {
            set_cache
                .entry(existing)
                .or_insert_with(|| nn_map[existing].iter().copied().collect());

            let candidate_set = &set_cache[&candidate];
            let existing_set = &set_cache[&existing];

            let shared = k2 - candidate_set.union(existing_set).count();
            max_overlap = max_overlap.max(shared);
        }

        if max_overlap <= max_shared {
            chosen.push(candidate);
        }
    }

    chosen
        .iter()
        .map(|&center| nn_map[center].as_slice())
        .collect()
}

pub fn aggregate_meta_cells(
    reader: &ParallelSparseReader,
    metacells: &[Vec<usize>],
    target_size: f32,
    n_genes: usize,
) -> CompressedSparseData<u32, f32> {
    let n_metacells = metacells.len();
    let mut all_data: Vec<u32> = Vec::new();
    let mut all_data_norm: Vec<f32> = Vec::new();
    let mut all_indices: Vec<usize> = Vec::new();
    let mut all_indptr: Vec<usize> = vec![0];

    const CHUNK_SIZE: usize = 1000;

    for chunk_start in (0..n_metacells).step_by(CHUNK_SIZE) {
        let chunk_end = (chunk_start + CHUNK_SIZE).min(n_metacells);
        let chunk = &metacells[chunk_start..chunk_end];

        let results: Vec<(Vec<usize>, Vec<u32>, Vec<f32>)> = chunk
            .par_iter()
            .map(|cell_idx| {
                let cells = reader.read_cells_parallel(&cell_idx);
                let mut gene_counts: FxHashMap<usize, u32> = FxHashMap::default();
                let mut library_size: u32 = 0;

                for cell in &cells {
                    for (idx, &count) in cell.indices.iter().zip(cell.data_raw.iter()) {
                        *gene_counts.entry(*idx as usize).or_insert(0) += count as u32;
                        library_size += count as u32;
                    }
                }

                // sort for CSR format
                let mut entries: Vec<(usize, u32)> = gene_counts.into_iter().collect();
                entries.sort_by_key(|(idx, _)| *idx);

                let indices: Vec<usize> = entries.iter().map(|(idx, _)| *idx).collect();
                let raw_counts: Vec<u32> = entries.iter().map(|(_, count)| *count).collect();
                let norm_counts: Vec<f32> = entries
                    .iter()
                    .map(|(_, count)| {
                        let norm = (*count as f32 / library_size as f32) * target_size;
                        (norm + 1.0).ln()
                    })
                    .collect();

                (indices, raw_counts, norm_counts)
            })
            .collect();

        for (indices, raw_counts, norm_counts) in results {
            all_indices.extend(indices);
            all_data.extend(raw_counts);
            all_data_norm.extend(norm_counts);
            all_indptr.push(all_indices.len());
        }
    }

    CompressedSparseData::new_csr(
        &all_data,
        &all_indices,
        &all_indptr,
        Some(&all_data_norm),
        (n_metacells, n_genes),
    )
}

////////////////////
// Main functions //
////////////////////

/// Pseudo-bulk data across cells based on cell indices
///
/// ### Params
///
/// * `f_path` - File path to the cell-based binary file.
/// * `cell_indices` - Slice of indices to pseudo-bulk.
/// * `verbose` - Controls verbosity of the function.
///
/// ### Return
///
/// Matrix of samples x genes pseudo-bulked.
pub fn get_pseudo_bulked_counts(
    f_path: &str,
    cell_indices: &[Vec<usize>],
    verbose: bool,
) -> Mat<f64> {
    let reader = ParallelSparseReader::new(f_path).unwrap();

    let n_genes = reader.get_header().total_genes;
    let n_groups = cell_indices.len();

    // Initialize matrix: groups (rows) × genes (columns)
    let mut result = Mat::zeros(n_groups, n_genes);

    for (group_idx, indices) in cell_indices.iter().enumerate() {
        let start_group = Instant::now();
        let chunks = reader.read_cells_parallel(indices);

        for chunk in chunks {
            for (value, &gene_idx) in chunk.data_raw.iter().zip(chunk.indices.iter()) {
                result[(group_idx, gene_idx as usize)] += *value as f64;
            }
        }

        if verbose && (group_idx + 1) % 10 == 0 {
            let elapsed = start_group.elapsed();
            let pct_complete = ((group_idx + 1) as f32 / n_groups as f32) * 100.0;
            println!(
                "Processed group {} out of {} (took {:.2?}, completed {:.1}%)",
                group_idx + 1,
                n_groups,
                elapsed,
                pct_complete
            );
        }
    }

    result
}
