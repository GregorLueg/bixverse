// use statrs::distribution::{Discrete, Hypergeometric};
use rustc_hash::FxHashSet;
use statrs::function::gamma::ln_gamma;

use crate::utils::utils_stats::calc_fdr;

///////////
// Types //
///////////

/// A type alias that can be returned by the par_iter() functions.
///
/// ### Fields
///
/// * `pval` - P-value
/// * `fdr` - FDR
/// * `odds_ratio` - Oddsratio
/// * `hits` - Number of success
/// * `gs_length` - Length of the gene set
pub struct HypergeomResult {
    pub pval: Vec<f64>,
    pub fdr: Vec<f64>,
    pub odds_ratio: Vec<f64>,
    pub hits: Vec<usize>,
    pub gs_length: Vec<usize>,
}

///////////////
// Functions //
///////////////

/// Calculate the p-value of a hypergeometric test.
///
/// ### Params
///
/// * `q` - Number of white balls drawn
/// * `m` - Number of white balls in the urn
/// * `n` - Number of black balls in the urn
/// * `k` - Number of balls drawn from the urn
///
/// ### Return
///
/// The p-value of the hypergeometric test
#[inline]
pub fn hypergeom_pval(q: usize, m: usize, n: usize, k: usize) -> f64 {
    if q == 0 {
        1.0
    } else {
        let population = m + n;

        // Always use logarithmic calculation to avoid numerical issues
        // Convert to f64 once at the start
        let (n_f, m_f, k_f) = (n as f64, m as f64, k as f64);
        let population_f = population as f64;

        // Calculate P(X > q) in log space
        let upper = k.min(m);

        // Use log space to compute probabilities for each value i > q
        let mut log_probs = Vec::new();
        for i in (q + 1)..=upper {
            let i_f = i as f64;

            // Calculate log(PMF(i)) using logarithms
            let log_pmf = ln_gamma(m_f + 1.0) - ln_gamma(i_f + 1.0) - ln_gamma(m_f - i_f + 1.0)
                + ln_gamma(n_f + 1.0)
                - ln_gamma(k_f - i_f + 1.0)
                - ln_gamma(n_f - (k_f - i_f) + 1.0)
                - (ln_gamma(population_f + 1.0)
                    - ln_gamma(k_f + 1.0)
                    - ln_gamma(population_f - k_f + 1.0));

            log_probs.push(log_pmf);
        }

        // If there are no probabilities to sum, return 0
        if log_probs.is_empty() {
            return 0.0;
        }

        // Use log-sum-exp trick to calculate sum without overflow
        let max_log_prob = log_probs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

        // Sum with adjustment to avoid numerical issues
        let mut sum = 0.0;
        for log_p in log_probs {
            sum += (log_p - max_log_prob).exp();
        }

        // Final result
        sum * max_log_prob.exp()
    }
}

/// Calculate odds ratios
///
/// ### Params
///
/// * `a1_b1` - In both gene set and target set
/// * `a0_b1` - In gene set, but not in target set
/// * `a1_b0` - In target set, but not in gene set
/// * `a0_b0` - Not in either
///
/// ### Return
///
/// The odds ratio. Pending values, can become infinity.
#[inline]
pub fn hypergeom_odds_ratio(a1_b1: usize, a0_b1: usize, a1_b0: usize, a0_b0: usize) -> f64 {
    (a1_b1 as f64 / a0_b1 as f64) / (a1_b0 as f64 / a0_b0 as f64)
}

/// Count the number of hits for the hypergeometric tests
///
/// ### Params
///
/// * `gene_set_list` - A slice of String vectors, representing the gene sets
///    you want to count the number of hits against
/// * `target_genes` - A string slice representing the target genes
///
/// ### Returns
///
/// A vector of hits, i.e., intersecting genes.
#[inline]
pub fn count_hits(
    gene_set_list: &[FxHashSet<String>],
    target_genes: &FxHashSet<String>,
) -> Vec<usize> {
    let hits: Vec<usize> = gene_set_list
        .iter()
        .map(|targets| {
            let intersection = targets.intersection(target_genes).count();
            intersection
        })
        .collect();

    hits
}

/// Helper function for the hypergeometric test
///
/// ### Params
///
/// - `target_genes` - The target genes for the test
/// - `gene_sets` - The list of vectors with the gene set genes
/// - `gene_universe` - Vector with the all genes of the universe
///
/// ### Returns
///
/// `HypergeomResult` - A tuple with the results.
pub fn hypergeom_helper(
    target_genes: &FxHashSet<String>,
    gene_sets: &[FxHashSet<String>],
    gene_universe: &[String],
) -> HypergeomResult {
    let gene_universe_length = gene_universe.len();
    let trials = target_genes.len();
    let gene_set_lengths = gene_sets.iter().map(|s| s.len()).collect::<Vec<usize>>();

    let hits = count_hits(gene_sets, target_genes);

    let pvals: Vec<f64> = hits
        .iter()
        .zip(gene_set_lengths.iter())
        .map(|(hit, gene_set_length)| {
            let q = *hit as i64 - 1;
            if q > 0 {
                hypergeom_pval(
                    q as usize,
                    *gene_set_length,
                    gene_universe_length - *gene_set_length,
                    trials,
                )
            } else {
                1.0
            }
        })
        .collect();

    let odds_ratios: Vec<f64> = hits
        .iter()
        .zip(gene_set_lengths.iter())
        .map(|(hit, gene_set_length)| {
            hypergeom_odds_ratio(
                *hit,
                *gene_set_length - *hit,
                trials - *hit,
                gene_universe_length - *gene_set_length - trials + *hit,
            )
        })
        .collect();

    let fdr = calc_fdr(&pvals);

    HypergeomResult {
        pval: pvals,
        fdr,
        odds_ratio: odds_ratios,
        hits,
        gs_length: gene_set_lengths,
    }
}

/// Helper function for the hypergeometric test
///
/// ### Params
///
/// - `res` - The `HypergeomResult` to filter.
/// - `min_overlap` - Optional minimum overlap in terms of hits.
/// - `fdr_threshold` - Optional threshold on the fdr.
///
/// ### Returns
///
/// `HypergeomResult` - The filtered results
pub fn filter_gse_results(
    res: HypergeomResult,
    min_overlap: Option<usize>,
    fdr_threshold: Option<f64>,
) -> (HypergeomResult, Vec<usize>) {
    let to_keep: Vec<usize> = (0..res.pval.len())
        .filter(|i| {
            if let Some(min_overlap) = min_overlap {
                if res.hits[*i] < min_overlap {
                    return false;
                }
            }
            if let Some(fdr_threshold) = fdr_threshold {
                if res.fdr[*i] > fdr_threshold {
                    return false;
                }
            }
            true
        })
        .collect();

    (
        HypergeomResult {
            pval: to_keep.iter().map(|i| res.pval[*i]).collect(),
            fdr: to_keep.iter().map(|i| res.fdr[*i]).collect(),
            odds_ratio: to_keep.iter().map(|i| res.odds_ratio[*i]).collect(),
            hits: to_keep.iter().map(|i| res.hits[*i]).collect(),
            gs_length: to_keep.iter().map(|i| res.gs_length[*i]).collect(),
        },
        to_keep.iter().map(|x| *x + 1).collect(),
    )
}
