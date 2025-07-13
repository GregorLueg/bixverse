// use statrs::distribution::{Discrete, Hypergeometric};
use rustc_hash::FxHashSet;
use statrs::function::gamma::ln_gamma;

///////////
// Types //
///////////

/// A type alias that can be returned by the par_iter() functions.
///
/// ### Fields
/// * `0` - P-value
/// * `1` - Odds ratio  
/// * `2` - Success counts
/// * `3` - Gene set size
pub type HypergeomResult = (Vec<f64>, Vec<f64>, Vec<usize>, Vec<usize>);

///////////////
// Functions //
///////////////

/// Calculate the p-value of a hypergeometric test.
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
pub fn hypergeom_odds_ratio(a1_b1: usize, a0_b1: usize, a1_b0: usize, a0_b0: usize) -> f64 {
    (a1_b1 as f64 / a0_b1 as f64) / (a1_b0 as f64 / a0_b0 as f64)
}

/// Count the number of hits for the hypergeometric tests
pub fn count_hits(gene_set_list: &[Vec<String>], target_genes: &[String]) -> Vec<usize> {
    let target_genes_hash: FxHashSet<_> = target_genes.iter().collect();
    let hits: Vec<usize> = gene_set_list
        .iter()
        .map(|s| {
            let s_hash: FxHashSet<_> = s.iter().collect();
            let intersection = s_hash.intersection(&target_genes_hash).count();
            intersection
        })
        .collect();

    hits
}

// /// Count the number of hits for the hypergeometric tests (against HashSets)
// pub fn count_hits_hash(gene_set_list: Vec<&HashSet<String>>, target_genes: &[String]) -> Vec<u64> {
//     let target_genes_hash: HashSet<String> = target_genes.iter().cloned().collect();
//     let hits: Vec<u64> = gene_set_list
//         .into_iter()
//         .map(|s| {
//             let intersection = s.intersection(&target_genes_hash).count() as u64;
//             intersection
//         })
//         .collect();
//     hits
// }

/// Helper function for the hypergeometric test
///
/// ### Parameters
///
/// - `target_genes` - The target genes for the test
/// - `gene_sets` - The list of vectors with the gene set genes
/// - `gene_universe` - Vector with the all genes of the universe
///
/// ### Returns
///
/// `HypergeomResult` - A tuple with the results.
pub fn hypergeom_helper(
    target_genes: &[String],
    gene_sets: &[Vec<String>],
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

    (pvals, odds_ratios, hits, gene_set_lengths)
}
