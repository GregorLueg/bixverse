use num_traits::Float;
use rayon::prelude::*;
use statrs::distribution::{Continuous, ContinuousCDF, Normal};
use statrs::function::gamma::ln_gamma;
use std::ops::{Add, Div};

///////////
// Types //
///////////

/////////////////////
// Enums | Helpers //
/////////////////////

#[derive(Clone, Debug)]
pub enum TestAlternative {
    /// Two sided test for the Z-score
    TwoSided,
    /// One-sided test for greater than
    Greater,
    /// One-sided test for lesser than
    Less,
}

/// Helper function to get the test alternative
///
/// ### Params
///
/// * `s` - String, type of test to run.
///
/// ### Returns
///
/// Option of the `TestAlternative`
pub fn get_test_alternative(s: &str) -> Option<TestAlternative> {
    match s.to_lowercase().as_str() {
        "greater" => Some(TestAlternative::Greater),
        "less" => Some(TestAlternative::Less),
        "twosided" => Some(TestAlternative::TwoSided),
        _ => None,
    }
}

///////////////
// Functions //
///////////////

/// Transform Z-scores into p-values (assuming normality).
///
/// ### Params
///
/// * `z_scores` - The Z scores to transform to p-values
///
/// ### Returns
///
/// The p-value vector based on the Z scores (two sided)
pub fn z_scores_to_pval(z_scores: &[f64], test_alternative: &str) -> Result<Vec<f64>, String> {
    let test_alternative = get_test_alternative(test_alternative)
        .ok_or_else(|| format!("Invalid Test alternative: {}", test_alternative))?;
    let normal = Normal::new(0.0, 1.0).unwrap();
    let res = z_scores
        .iter()
        .map(|&z| match test_alternative {
            TestAlternative::TwoSided => {
                let abs_z = z.abs();
                if abs_z > 6.0 {
                    let pdf = normal.pdf(abs_z);
                    let p = pdf / abs_z * (1.0 - 1.0 / (abs_z * abs_z));
                    2.0 * p
                } else {
                    2.0 * (1.0 - normal.cdf(abs_z))
                }
            }
            TestAlternative::Greater => {
                if z > 6.0 {
                    let pdf = normal.pdf(z);
                    pdf / z * (1.0 - 1.0 / (z * z))
                } else {
                    1.0 - normal.cdf(z)
                }
            }
            TestAlternative::Less => {
                if z < -6.0 {
                    let abs_z = z.abs();
                    let pdf = normal.pdf(abs_z);
                    pdf / abs_z * (1.0 - 1.0 / (abs_z * abs_z))
                } else {
                    normal.cdf(z)
                }
            }
        })
        .collect();
    Ok(res)
}

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

/// Calculate the FDR
///
/// ### Params
///
/// * `pvals` - P-values for which to calculate the FDR
///
/// ### Returns
///
/// The calculated FDRs
pub fn calc_fdr(pvals: &[f64]) -> Vec<f64> {
    let n = pvals.len();
    let n_f64 = n as f64;

    let mut indexed_pval: Vec<(usize, f64)> =
        pvals.par_iter().enumerate().map(|(i, &x)| (i, x)).collect();

    indexed_pval
        .sort_unstable_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));

    let adj_pvals_tmp: Vec<f64> = indexed_pval
        .par_iter()
        .enumerate()
        .map(|(i, (_, p))| (n_f64 / (i + 1) as f64) * p)
        .collect();

    let mut current_min = adj_pvals_tmp[n - 1].min(1.0);
    let mut monotonic_adj = vec![current_min; n];

    for i in (0..n - 1).rev() {
        current_min = current_min.min(adj_pvals_tmp[i]).min(1.0);
        monotonic_adj[i] = current_min;
    }

    let mut adj_pvals = vec![0.0; n];

    for (i, &(original_idx, _)) in indexed_pval.iter().enumerate() {
        adj_pvals[original_idx] = monotonic_adj[i];
    }

    adj_pvals
}

/// Get the median
///
/// ### Params
///
/// * `x` - The slice for which to calculate the median for.
///
/// ### Results
///
/// The median (if the vector is not empty)
pub fn median<T>(x: &[T]) -> Option<T>
where
    T: Clone + PartialOrd + Add<Output = T> + Div<T, Output = T> + From<u8>,
{
    if x.is_empty() {
        return None;
    }

    let mut data = x.to_vec();
    let len = data.len();

    if len.is_multiple_of(2) {
        let (_, median1, right) =
            data.select_nth_unstable_by(len / 2 - 1, |a, b| a.partial_cmp(b).unwrap());
        let median2 = right
            .iter()
            .min_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap();
        Some((median1.clone() + median2.clone()) / T::from(2))
    } else {
        let (_, median, _) = data.select_nth_unstable_by(len / 2, |a, b| a.partial_cmp(b).unwrap());
        Some(median.clone())
    }
}

/// Calculate the MAD
///
/// ### Params
///
/// * `x` - Slice for which to calculate the MAD for
///
/// ### Results
///
/// The MAD of the slice.
pub fn mad(x: &[f64]) -> Option<f64> {
    if x.is_empty() {
        return None;
    }

    let median_val = median(x)?; // Early return if median is None
    let deviations: Vec<f64> = x.iter().map(|&x| (x - median_val).abs()).collect();

    median(&deviations)
}

///////////
// Other //
///////////

/// Logit function
///
/// ### Params
///
/// * `p` - Probability value (must be in (0, 1))
///
/// ### Returns
///
/// Log-odds: ln(p / (1-p))
pub fn logit<F: Float>(p: F) -> F {
    (p / (F::one() - p)).ln()
}

/// Inverse logit (sigmoid) function
///
/// ### Params
///
/// * `q` - Log-odds value
///
/// ### Returns
///
/// Probability: exp(q) / (1 + exp(q))
pub fn inv_logit<F: Float>(q: F) -> F {
    q.exp() / (F::one() + q.exp())
}
