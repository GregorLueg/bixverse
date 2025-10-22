use extendr_api::prelude::*;
use rand::prelude::*;
use rayon::prelude::*;

use crate::core::base::stats::*;
use crate::core::base::utils::{col_means, col_sds};
use crate::utils::general::split_vector_randomly;
use crate::utils::r_rust_interface::*;

/// Fast AUC calculation
///
/// @description This function calculates rapidly AUCs based on an approximation.
///
/// @param pos_scores The scores of your hits.
/// @param neg_scores The scores of your non-hits.
/// @param iters Number of iterations to run the function for.
/// Recommended size: 10000L.
/// @param seed Seed.
///
/// @return The AUC.
///
/// @export
#[extendr]
fn rs_fast_auc(pos_scores: &[f64], neg_scores: &[f64], iters: usize, seed: u64) -> f64 {
    let mut rng = StdRng::seed_from_u64(seed);
    let mut count = 0;

    for _ in 0..iters {
        let pos_sample = pos_scores.choose(&mut rng).unwrap();
        let neg_sample = neg_scores.choose(&mut rng).unwrap();
        if pos_sample > neg_sample {
            count += 1;
        }
    }

    count as f64 / iters as f64
}

/// Create random AUCs
///
/// @description This function creates a random set of AUCs based on a score
/// vector and a size of the positive set. This can be used for permutation-
/// based estimation of Z-scores and subsequently p-values.
///
/// @param score_vec The overall vector of scores.
/// @param size_pos The size of the hits represented in the score_vec.
/// @param random_iters Number of random AUCs to generate.
/// @param auc_iters Number of random iterations to approximate the AUCs.
/// Recommended size: 10000L.
/// @param seed Seed.
///
/// @return A vector of random AUCs based the score vector and size of the
/// positive set.
///
/// @export
#[extendr]
fn rs_create_random_aucs(
    score_vec: &[f64],
    size_pos: usize,
    random_iters: usize,
    auc_iters: usize,
    seed: u64,
) -> Vec<f64> {
    let iter_vec: Vec<usize> = (0..random_iters).collect();

    let random_aucs: Vec<_> = iter_vec
        .par_iter()
        .map(|x| {
            let scores = split_vector_randomly(score_vec, size_pos, *x as u64 + seed);

            rs_fast_auc(&scores.0, &scores.1, auc_iters, *x as u64 + 1 + seed)
        })
        .collect();

    random_aucs
}

/// Calculate the Hedge's G effect
///
/// @description Calculates the Hedge's G effect for two sets of matrices. The
/// function assumes that rows = samples and columns = features.
/// WARNING! Incorrect use can cause kernel crashes. Wrapper around the Rust
/// functions with type checks are provided in the package.
///
/// @param mat_a The matrix of samples and features in grp A for which to
/// calculate the Hedge's G effect.
/// @param mat_b The matrix of samples and features in grp B for which to
/// calculate the Hedge's G effect.
/// @param small_sample_correction Shall the small sample correction be applied.
///
/// @return Returns the harmonic sum according to the OT calculation.
///
/// @export
#[extendr]
fn rs_hedges_g(mat_a: RMatrix<f64>, mat_b: RMatrix<f64>, small_sample_correction: bool) -> List {
    let mat_a = r_matrix_to_faer(&mat_a);
    let mat_b = r_matrix_to_faer(&mat_b);

    let n_a = mat_a.nrows();
    let n_b = mat_b.nrows();

    let mean_a = col_means(mat_a);
    let mean_b = col_means(mat_b);

    let std_a = col_sds(mat_a);
    let std_b = col_sds(mat_b);

    let (es, se): EffectSizeRes = hedge_g_effect(
        &mean_a,
        &mean_b,
        &std_a,
        &std_b,
        n_a,
        n_b,
        small_sample_correction,
    );

    list!(effect_sizes = es, standard_errors = se)
}

/// Calculate a BH-based FDR
///
/// @description Rust implementation that will be faster if you have an
/// terrifying amount of p-values to adjust.
///
/// @param pvals Numeric vector. The p-values you wish to adjust.
///
/// @return The Benjamini-Hochberg adjusted p-values.
///
/// @export
#[extendr]
fn rs_fdr_adjustment(pvals: &[f64]) -> Vec<f64> {
    calc_fdr(pvals)
}

/// Calculate the hypergeometric rest in Rust
///
/// @param q Number of white balls drawn out of urn.
/// @param m Number of white balls in the urn.
/// @param n Number of black balls in the urn.
/// @param k The number of balls drawn out of the urn.
///
/// @return P-value (with lower.tail set to False)
///
/// @export
#[extendr]
fn rs_phyper(q: usize, m: usize, n: usize, k: usize) -> f64 {
    hypergeom_pval(q, m, n, k)
}

/// Calculate MAD outlier detection in Rust.
///
/// @param x Numerical vector to test.
/// @param threshold Numeric. Number of MADs in either direction that is
/// acceptable.
/// @param direction String. One of `c("below", "above", "twosided")`. Shall
/// the outlier direction be done for values below the threshold, above the
/// threshold or in both directions.
///
/// @return A list with the following items:
/// \itemize{
///  \item outlier - Boolean vector if element is an outlier
///  \item threshold - Applied final threshold
/// }
///
/// @details
/// Should you provide too short vectors, the function will return an empty
/// boolean and a threshold of 0.
///
/// @export
#[extendr]
fn rs_mad_outlier(x: &[f64], threshold: f64, direction: &str) -> extendr_api::Result<List> {
    let direction = get_outlier_type(direction)
        .ok_or_else(|| format!("Invalid direction type: {}", direction))?;

    let (outliers, threshold) = mad_outlier(x, threshold, direction);

    Ok(list!(outlier = outliers, threshold = threshold))
}

extendr_module! {
    mod r_stats;
    fn rs_fast_auc;
    fn rs_create_random_aucs;
    fn rs_hedges_g;
    fn rs_fdr_adjustment;
    fn rs_phyper;
    fn rs_mad_outlier;
}
