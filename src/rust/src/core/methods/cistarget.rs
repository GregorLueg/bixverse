use extendr_api::{Attributes, List, Robj};
use faer::MatRef;
use rayon::prelude::*;

use crate::utils::general::standard_deviation;

////////////////
// Structures //
////////////////

/// Results structure for a single enriched motif
///
/// ### Fields
///
/// * `motif_idx` - Index of the motif
/// * `nes` - Normalised enrichment score
/// * `auc` - Area under the curve
/// * `rank_at_max` - Rank position at maximum enrichment
/// * `n_enriched` - Number of enriched genes
/// * `enriched_gene_indices` - Indices of enriched genes in the leading edge
#[derive(Debug, Clone)]
pub struct MotifEnrichment {
    pub motif_idx: usize,
    pub nes: f64,
    pub auc: f64,
    pub rank_at_max: u32,
    pub n_enriched: usize,
    pub enriched_gene_indices: Vec<usize>,
}

///////////
// Enums //
///////////

/// Enum for the ICA types
#[derive(Clone, Debug)]
pub enum RccType {
    /// Use the `ICisTarget` type for the RCC calculation
    ICisTarget,
    /// Use the `Approx` type for the RCC calculation
    Approx,
}

/// Parse the RCC calculation type
///
/// ### Params
///
/// * `s` - The string to parse
///
/// ### Returns
///
/// Option of the `RccType`
pub fn parse_rcc_type(s: &str) -> Option<RccType> {
    match s.to_lowercase().as_str() {
        "icistarget" => Some(RccType::ICisTarget),
        "approx" => Some(RccType::Approx),
        _ => None,
    }
}

/////////////////////////
// Helpers - AUC & NES //
/////////////////////////

/// Calculate AUC for a single motif's gene ranks
///
/// ### Params
///
/// * `gene_ranks` - Slice of gene ranks for this motif
/// * `auc_threshold` - Maximum rank to consider for AUC calculation
///
/// ### Returns
///
/// Normalised AUC score
fn calculate_auc_single(gene_ranks: &[i32], auc_threshold: i32) -> f64 {
    let mut filtered = gene_ranks
        .iter()
        .copied()
        .filter(|&r| r < auc_threshold)
        .collect::<Vec<i32>>();

    if filtered.is_empty() {
        return 0.0;
    }

    filtered.sort_unstable();
    filtered.push(auc_threshold);

    let n_genes = gene_ranks.len() as f64;
    let max_auc = auc_threshold as f64 * n_genes;

    let auc = filtered
        .windows(2)
        .enumerate()
        .map(|(i, w)| (w[1] - w[0]) as f64 * (i + 1) as f64)
        .sum::<f64>();

    auc / max_auc
}

/// Calculate AUC scores for all motifs given a gene set
///
/// ### Params
///
/// * `matrix` - Matrix reference containing gene-motif ranks
/// * `gene_indices` - Indices of genes in the gene set
/// * `auc_threshold` - Maximum rank to consider for AUC calculation
///
/// ### Returns
///
/// Vector of AUC scores for each motif
pub fn calculate_aucs(
    matrix: MatRef<'_, i32>,
    gene_indices: &[usize],
    auc_threshold: i32,
) -> Vec<f64> {
    let n_motifs = matrix.ncols();

    (0..n_motifs)
        .into_par_iter()
        .map(|motif_idx| {
            let gene_ranks: Vec<i32> = gene_indices
                .iter()
                .map(|&row_idx| *matrix.get(row_idx, motif_idx))
                .collect();

            calculate_auc_single(&gene_ranks, auc_threshold)
        })
        .collect()
}

/// Calculate Normalised Enrichment Scores from AUC values
///
/// ### Params
///
/// * `auc_scores` - Slice of AUC scores
///
/// ### Returns
///
/// Vector of z-score normalised enrichment scores
pub fn calculate_nes(auc_scores: &[f64]) -> Vec<f64> {
    let n = auc_scores.len() as f64;
    let mean = auc_scores.iter().sum::<f64>() / n;

    let std = standard_deviation(auc_scores);

    if std == 0.0 {
        return vec![0.0; auc_scores.len()];
    }

    auc_scores
        .par_iter()
        .map(|&auc| (auc - mean) / std)
        .collect()
}

/////////
// RCC //
/////////

/// Calculate rolling mean with centre alignment and extend fill
///
/// ### Params
///
/// * `values` - Slice of values to smooth
/// * `window` - Window size for rolling mean
///
/// ### Returns
///
/// Vector of smoothed values
fn rolling_mean(values: &[f64], window: usize) -> Vec<f64> {
    let n = values.len();
    let mut result = Vec::with_capacity(n);
    let half_window = window / 2;
    let mut sum = 0.0;
    let mut count = 0;

    // Initial window
    let start = 0;
    let end = half_window.min(n);
    for &v in &values[start..end] {
        sum += v;
        count += 1;
    }
    result.push(sum / count as f64);

    // Slide window
    for i in 1..n {
        let new_start = i.saturating_sub(half_window);
        let new_end = (i + half_window + 1).min(n);
        let old_start = (i - 1).saturating_sub(half_window);
        let old_end = (i - 1 + half_window + 1).min(n);

        // Remove values that left the window
        for &v in &values[old_start..new_start] {
            sum -= v;
            count -= 1;
        }

        // Add values entering the window
        for &v in &values[old_end..new_end] {
            sum += v;
            count += 1;
        }

        result.push(sum / count as f64);
    }

    result
}

/// Calculate Recovery Curve for a single motif
///
/// ### Params
///
/// * `gene_ranks` - Slice of gene ranks for this motif
/// * `max_rank` - Maximum rank to consider
///
/// ### Returns
///
/// Vector where position i equals the number of genes recovered by rank i
fn calculate_rcc(gene_ranks: &[i32], max_rank: i32) -> Vec<usize> {
    let mut rcc = vec![0; max_rank as usize];
    let mut filtered: Vec<i32> = gene_ranks
        .iter()
        .copied()
        .filter(|&r| r > 0 && r < max_rank)
        .collect();

    if filtered.is_empty() {
        return rcc;
    }

    filtered.sort_unstable();
    filtered.push(max_rank);

    // Build step function: at each rank, count cumulative genes
    let mut prev_rank = 0;
    for (gene_count, &rank) in filtered.iter().enumerate() {
        // Fill from prev_rank to current rank with gene_count
        #[allow(clippy::needless_range_loop)]
        for pos in prev_rank..(rank as usize).min(max_rank as usize) {
            rcc[pos] = gene_count;
        }
        prev_rank = rank as usize;
    }

    rcc
}

/// Calculate RCC for all motifs
///
/// ### Params
///
/// * `matrix` - Matrix reference containing gene-motif ranks
/// * `gene_indices` - Indices of genes in the gene set
/// * `max_rank` - Maximum rank to consider
///
/// ### Returns
///
/// Vector of RCC curves for each motif
fn calculate_all_rccs(
    matrix: MatRef<'_, i32>,
    gene_indices: &[usize],
    max_rank: i32,
) -> Vec<Vec<usize>> {
    let n_motifs = matrix.ncols();

    (0..n_motifs)
        .into_par_iter()
        .map(|motif_idx| {
            let gene_ranks: Vec<i32> = gene_indices
                .iter()
                .map(|&idx| *matrix.get(idx, motif_idx))
                .collect();

            calculate_rcc(&gene_ranks, max_rank)
        })
        .collect()
}

/// Calculate mean+2*sd thresholds from pre-calculated RCCs (exact method)
///
/// ### Params
///
/// * `all_rccs` - Pre-calculated RCC curves for all motifs
/// * `max_rank` - Maximum rank to consider
///
/// ### Returns
///
/// Vector of threshold values at each rank position
fn calculate_rcc_thresholds(all_rccs: &[Vec<usize>], max_rank: i32) -> Vec<f64> {
    let max_rank_usize = max_rank as usize;
    let n_motifs = all_rccs.len() as f64;
    let mut rcc_m2sd = vec![0.0; max_rank_usize];

    for rank_pos in 0..max_rank_usize {
        let mut sum = 0.0;
        let mut sum_sq = 0.0;

        for rcc in all_rccs {
            let val = rcc[rank_pos] as f64;
            sum += val;
            sum_sq += val * val;
        }

        let mean = sum / n_motifs;
        let variance = (sum_sq / n_motifs) - (mean * mean);
        let sd = variance.sqrt();

        rcc_m2sd[rank_pos] = mean + 2.0 * sd;
    }

    rcc_m2sd
}

/// Calculate mean+2*sd thresholds using approximate method
///
/// ### Params
///
/// * `matrix` - Matrix reference containing gene-motif ranks
/// * `gene_indices` - Indices of genes in the gene set
/// * `max_rank` - Maximum rank to consider
/// * `n_mean` - Window size for smoothing
///
/// ### Returns
///
/// Vector of smoothed threshold values at each rank position
fn calculate_rcc_thresholds_aprox(
    matrix: MatRef<'_, i32>,
    gene_indices: &[usize],
    max_rank: i32,
    n_mean: usize,
) -> Vec<f64> {
    let n_motifs = matrix.ncols();
    let n_genes = gene_indices.len();
    let max_rank_extra = (max_rank as usize) + n_mean;

    let mut global_mat = vec![vec![0usize; max_rank_extra]; n_genes + 1];

    for m_idx in 0..n_motifs {
        let mut ranks: Vec<i32> = gene_indices
            .iter()
            .map(|&idx| *matrix.get(idx, m_idx))
            .filter(|&r| r > 0 && r <= max_rank_extra as i32)
            .collect();

        if ranks.is_empty() {
            continue;
        }

        ranks.sort_unstable();

        for (gene_count, &rank) in ranks.iter().enumerate() {
            if (rank as usize) < max_rank_extra {
                global_mat[gene_count + 1][rank as usize] += 1;
            }
        }
    }

    let mut rcc_mean_raw = vec![0.0; max_rank_extra];
    let mut rcc_sd_raw = vec![0.0; max_rank_extra];

    for rank_pos in 0..max_rank_extra {
        let total_count: usize = global_mat.iter().map(|row| row[rank_pos]).sum();

        if total_count == 0 {
            if rank_pos > 0 {
                rcc_mean_raw[rank_pos] = rcc_mean_raw[rank_pos - 1];
                rcc_sd_raw[rank_pos] = rcc_sd_raw[rank_pos - 1];
            }
            continue;
        }

        let total_count_f64 = total_count as f64;

        let mut sum = 0.0;
        let mut sum_sq = 0.0;

        for (num_genes, count) in global_mat.iter().enumerate() {
            let c = count[rank_pos] as f64;
            let ng = num_genes as f64;
            sum += ng * c;
            sum_sq += ng * ng * c;
        }

        let mean = sum / total_count_f64;
        let variance = (sum_sq / total_count_f64) - (mean * mean);

        rcc_mean_raw[rank_pos] = mean;
        rcc_sd_raw[rank_pos] = variance.sqrt();
    }

    let mut rcc_mean = vec![0.0; max_rank as usize];
    let mut rcc_sd = vec![0.0; max_rank as usize];

    let smoothed_mean = rolling_mean(&rcc_mean_raw, n_mean);
    let smoothed_sd = rolling_mean(&rcc_sd_raw, n_mean);

    for i in 0..(max_rank as usize) {
        if i < 5 {
            rcc_mean[i] = rcc_mean_raw[i];
            rcc_sd[i] = rcc_sd_raw[i];
        } else if i < smoothed_mean.len() {
            rcc_mean[i] = smoothed_mean[i];
            rcc_sd[i] = smoothed_sd[i];
        }
    }

    rcc_mean
        .iter()
        .zip(rcc_sd.iter())
        .map(|(m, s)| m + 2.0 * s)
        .collect()
}

/// Find the leading edge for a motif
///
/// ### Params
///
/// * `matrix` - Matrix reference containing gene-motif ranks
/// * `gene_indices` - Indices of genes in the gene set
/// * `motif_idx` - Index of the motif to analyse
/// * `rcc_m2sd` - Pre-calculated threshold values
/// * `all_rccs` - Optional pre-calculated RCC curves (for exact method)
///
/// ### Returns
///
/// Tuple of (rank at maximum enrichment, number of enriched genes, enriched gene indices)
fn find_leading_edge(
    matrix: MatRef<'_, i32>,
    gene_indices: &[usize],
    motif_idx: usize,
    rcc_m2sd: &[f64],
    all_rccs: Option<&Vec<Vec<usize>>>,
) -> (u32, usize, Vec<usize>) {
    let current_rcc = match all_rccs {
        Some(rccs) => &rccs[motif_idx],
        None => &calculate_rcc(
            &gene_indices
                .iter()
                .map(|&idx| *matrix.get(idx, motif_idx))
                .collect::<Vec<_>>(),
            rcc_m2sd.len() as i32,
        ),
    };

    let mut max_enrichment = f64::NEG_INFINITY;
    let mut best_rank = 0;

    for (rank_pos, &threshold) in rcc_m2sd.iter().enumerate() {
        if rank_pos < current_rcc.len() {
            let enrichment = current_rcc[rank_pos] as f64 - threshold;
            if enrichment > max_enrichment {
                max_enrichment = enrichment;
                best_rank = rank_pos;
            }
        }
    }

    let rank_at_max = best_rank as u32;

    let enriched: Vec<usize> = gene_indices
        .iter()
        .filter(|&&idx| {
            let rank = *matrix.get(idx, motif_idx);
            rank > 0 && rank <= rank_at_max as i32
        })
        .copied()
        .collect();

    let n_enriched = enriched.len();

    (rank_at_max, n_enriched, enriched)
}

//////////
// Main //
//////////

/// Convert motif enrichments to R list format
///
/// ### Params
///
/// * `enrichments` - Slice of motif enrichment results
///
/// ### Returns
///
/// R list containing `motif_idx`, `nes`, `auc`, `rank_at_max`, `n_enriched`,
/// and `leading_edge`
pub fn motif_enrichments_to_r_list(enrichments: &[MotifEnrichment]) -> List {
    let mut result = List::new(6);

    let motif_idx: Vec<i32> = enrichments
        .iter()
        .map(|m| (m.motif_idx + 1) as i32)
        .collect();
    let nes: Vec<f64> = enrichments.iter().map(|m| m.nes).collect();
    let auc: Vec<f64> = enrichments.iter().map(|m| m.auc).collect();
    let rank_at_max: Vec<i32> = enrichments.iter().map(|m| m.rank_at_max as i32).collect();
    let n_enriched: Vec<i32> = enrichments.iter().map(|m| m.n_enriched as i32).collect();

    let mut leading_edge = List::new(enrichments.len());
    for (j, motif) in enrichments.iter().enumerate() {
        let genes: Vec<i32> = motif
            .enriched_gene_indices
            .iter()
            .map(|&idx| (idx + 1) as i32)
            .collect();
        leading_edge.set_elt(j, Robj::from(genes)).unwrap();
    }

    result.set_elt(0, Robj::from(motif_idx)).unwrap();
    result.set_elt(1, Robj::from(nes)).unwrap();
    result.set_elt(2, Robj::from(auc)).unwrap();
    result.set_elt(3, Robj::from(rank_at_max)).unwrap();
    result.set_elt(4, Robj::from(n_enriched)).unwrap();
    result.set_elt(5, Robj::from(leading_edge)).unwrap();

    result
        .set_names(&[
            "motif_idx",
            "nes",
            "auc",
            "rank_at_max",
            "n_enriched",
            "leading_edge",
        ])
        .unwrap();

    result
}

/// Process a single gene set - calculate AUC, NES, and find significant motifs
///
/// ### Params
///
/// * `matrix` - Matrix reference containing gene-motif ranks
/// * `gene_indices` - Indices of genes in the gene set
/// * `auc_threshold` - Maximum rank to consider for AUC calculation
/// * `nes_threshold` - Minimum NES value for significance
/// * `max_rank` - Maximum rank to consider for RCC
/// * `method` - Method to use ("icistarget" for exact, other for approximate)
/// * `n_mean` - Window size for smoothing (used in approximate method)
///
/// ### Returns
///
/// Vector of significant motif enrichments
pub fn process_gene_set(
    matrix: MatRef<'_, i32>,
    gene_indices: &[usize],
    auc_threshold: i32,
    nes_threshold: f64,
    max_rank: i32,
    method: &RccType,
    n_mean: usize,
) -> Vec<MotifEnrichment> {
    let auc_scores = calculate_aucs(matrix, gene_indices, auc_threshold);

    let nes_scores = calculate_nes(&auc_scores);

    let (rcc_m2sd, all_rccs) = match method {
        RccType::ICisTarget => {
            let all_rccs = calculate_all_rccs(matrix, gene_indices, max_rank);
            let rcc_m2sd = calculate_rcc_thresholds(&all_rccs, max_rank);
            (rcc_m2sd, Some(all_rccs))
        }
        RccType::Approx => {
            let rcc_m2sd = calculate_rcc_thresholds_aprox(matrix, gene_indices, max_rank, n_mean);
            (rcc_m2sd, None)
        }
    };

    nes_scores
        .iter()
        .enumerate()
        .filter(|(_, &nes)| nes >= nes_threshold)
        .map(|(motif_idx, &nes)| {
            let (rank_at_max, n_enriched, enriched_gene_indices) = find_leading_edge(
                matrix,
                gene_indices,
                motif_idx,
                &rcc_m2sd,
                all_rccs.as_ref(),
            );

            MotifEnrichment {
                motif_idx,
                nes,
                auc: auc_scores[motif_idx],
                rank_at_max,
                n_enriched,
                enriched_gene_indices,
            }
        })
        .collect()
}
