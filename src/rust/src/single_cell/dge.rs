use rayon::prelude::*;
use std::time::Instant;

use crate::core::base::stats::{calc_fdr, z_scores_to_pval};
use crate::core::data::sparse_io::*;
use crate::single_cell::fast_ranking::rank_csr_chunk_vec;

////////////////
// Structures //
////////////////

/// Structure to store the Mann Whitney U-based DGE results in
pub struct DgeMannWhitneyRes {
    pub lfc: Vec<f32>,
    pub prop1: Vec<f32>,
    pub prop2: Vec<f32>,
    pub z_scores: Vec<f64>,
    pub p_vals: Vec<f64>,
    pub fdr: Vec<f64>,
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
        for (&gene_idx, &norm_val) in cell.col_indices.iter().zip(cell.data_norm.iter()) {
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

    (u1 - mean) / variance.sqrt()
}

////////////////////
// Main functions //
////////////////////

pub fn calculate_dge_grps_mann_whitney(
    f_path: &str,
    grp_1_indices: &[usize],
    grp_2_indices: &[usize],
    min_proportion: f32,
    verbose: bool,
) -> DgeMannWhitneyRes {
    let start_read = Instant::now();

    let reader = ParallelSparseReader::new(f_path).unwrap();
    let no_genes = reader.get_header().total_genes;

    let mut cell_chunks_1: Vec<CsrCellChunk> = reader.read_cells_parallel(grp_1_indices);
    let mut cell_chunks_2: Vec<CsrCellChunk> = reader.read_cells_parallel(grp_2_indices);

    let end_read = start_read.elapsed();

    if verbose {
        println!("Load in data: {:.2?}", end_read);
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

    let p_vals = z_scores_to_pval(&z_scores);
    let fdr = calc_fdr(&p_vals);

    let end_calculations = start_calculations.elapsed();

    if verbose {
        println!("Finished DGE calculations: {:.2?}", end_calculations);
    }

    DgeMannWhitneyRes {
        lfc: log_fc,
        prop1,
        prop2,
        z_scores,
        p_vals,
        fdr,
    }
}
