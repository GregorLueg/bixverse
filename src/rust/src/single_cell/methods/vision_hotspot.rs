use faer::MatRef;
use rand::rngs::SmallRng;
use rand::seq::SliceRandom;
use rand::SeedableRng;
use rayon::prelude::*;
use std::time::Instant;

use crate::core::base::utils::rank_matrix_col;

/// Calculate Geary's C for a single signature (pathway)
///
/// This implements a modified approach akin to VISION.
///
/// ### Params
///
/// * `scores` - Ranked signature scores for all cells
/// * `knn_indices` - KNN indices matrix (cells x k)
/// * `knn_weights` - KNN weights matrix (cells x k)
///
/// ### Returns
///
/// Geary's C statistic
fn geary_c(scores: &[f64], knn_indices: &[Vec<usize>], knn_weights: &[Vec<f32>]) -> f64 {
    let n = scores.len();

    let mean: f64 = scores.iter().sum::<f64>() / n as f64;
    let variance: f64 = scores.iter().map(|x| (x - mean).powi(2)).sum::<f64>();

    if variance == 0.0 {
        return 0.0;
    }

    let mut numerator = 0.0;
    let mut total_weight = 0.0;

    // unsafe unchecked access in hot loop
    for (i, (indices, weights)) in knn_indices.iter().zip(knn_weights.iter()).enumerate() {
        let xi = unsafe { *scores.get_unchecked(i) };
        for (&j, &w) in indices.iter().zip(weights.iter()) {
            let xj = unsafe { *scores.get_unchecked(j) };
            numerator += w as f64 * (xi - xj).powi(2);
            total_weight += w as f64;
        }
    }

    let norm = 2.0 * total_weight * variance / (n as f64 - 1.0);

    numerator / norm
}

/// Calculate KNN weights using exponential kernel
///
/// ### Params
///
/// * `knn_indices` - KNN indices (cells x k)
/// * `knn_distances` - KNN squared distances (cells x k) - use Euclidean here!
///
/// ### Returns
///
/// KNN weights matrix (cells x k)
fn calc_knn_weights(knn_indices: &[Vec<usize>], knn_distances: &[Vec<f32>]) -> Vec<Vec<f32>> {
    knn_indices
        .par_iter() // Parallel iterator
        .zip(knn_distances.par_iter())
        .map(|(indices, distances)| {
            if distances.is_empty() {
                return vec![];
            }

            let sigma_sq = distances.last().copied().unwrap_or(1.0);

            if sigma_sq == 0.0 {
                return vec![1.0; indices.len()];
            }

            distances
                .iter()
                .map(|&d_sq| (-d_sq / sigma_sq).exp())
                .collect()
        })
        .collect()
}

/// Calculate VISION local autocorrelation scores
///
/// ### Params
///
/// * `pathway_scores` - Matrix of VISION scores (cells x pathways)
/// * `knn_indices` - KNN indices from embedding (cells x k)
/// * `knn_distances` - KNN squared distances (cells x k)
/// * `n_perm` - Number of permutations for significance testing
/// * `verbose` - Print progress
///
/// ### Returns
///
/// Tuple of (consistency_scores, p_values) for each pathway
pub fn calc_local_autocorrelation(
    pathway_scores: MatRef<f64>,
    knn_indices: Vec<Vec<usize>>,
    knn_distances: Vec<Vec<f32>>,
    n_perm: usize,
    seed: usize,
    verbose: bool,
) -> (Vec<f64>, Vec<f64>) {
    let start = Instant::now();
    let ranked_scores = rank_matrix_col(&pathway_scores);
    let knn_weights = calc_knn_weights(&knn_indices, &knn_distances);

    if verbose {
        println!("Computed KNN weights: {:.2?}", start.elapsed());
    }

    let n_pathways = ranked_scores.ncols();
    let start_geary = Instant::now();

    // parallel processing
    let results: Vec<_> = (0..n_pathways)
        .into_par_iter()
        .map(|pathway_idx| {
            let scores: Vec<f64> = pathway_scores.col(pathway_idx).iter().copied().collect();
            let ranks: Vec<f64> = ranked_scores.col(pathway_idx).iter().copied().collect();
            let observed_c = geary_c(&ranks, &knn_indices, &knn_weights);
            let observed_consistency = 1.0 - observed_c;

            let mut rng = SmallRng::seed_from_u64((seed + pathway_idx) as u64);
            let mut num_greater = 0;

            let mut indices: Vec<usize> = (0..scores.len()).collect();

            for _ in 0..n_perm {
                indices.shuffle(&mut rng);
                let permuted_scores: Vec<f64> = indices.iter().map(|&i| scores[i]).collect();

                // Re-rank the permuted scores
                let mut paired: Vec<(f64, usize)> = permuted_scores
                    .iter()
                    .enumerate()
                    .map(|(i, &v)| (v, i))
                    .collect();
                paired.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
                let mut permuted_ranks = vec![0.0; permuted_scores.len()];
                for (rank, &(_val, idx)) in paired.iter().enumerate() {
                    permuted_ranks[idx] = rank as f64;
                }

                let perm_c = geary_c(&permuted_ranks, &knn_indices, &knn_weights);
                let perm_consistency = 1.0 - perm_c;

                if perm_consistency >= observed_consistency {
                    num_greater += 1;
                }
            }

            let p_value = (num_greater + 1) as f64 / (n_perm + 1) as f64;

            (observed_consistency, p_value)
        })
        .collect();

    if verbose {
        println!(
            "Calculated autocorrelation scores: {:.2?}",
            start_geary.elapsed()
        );
    }

    let (consistency_scores, p_values): (Vec<_>, Vec<_>) = results.into_iter().unzip();

    (consistency_scores, p_values)
}
