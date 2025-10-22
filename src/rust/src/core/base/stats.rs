use faer::{linalg::solvers::DenseSolveCore, Mat};
use rand::prelude::*;
use rayon::prelude::*;
use statrs::distribution::FisherSnedecor;
use statrs::distribution::{Continuous, ContinuousCDF, Normal};
use statrs::function::gamma::ln_gamma;
use std::ops::{Add, Div};

use crate::assert_same_len;

///////////
// Types //
///////////

/// A type alias representing effect size results
///
/// ### Fields
///
/// * `0` - The calculated effect sizes
/// * `1` - The corresponding standard errors
pub type EffectSizeRes = (Vec<f64>, Vec<f64>);

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

#[derive(Clone, Debug)]
pub enum OutlierDirection {
    /// Check if outlier is below the threshold
    Below,
    /// Check if outlier is above the threshold
    Above,
    /// Check if outlier is below OR above the thresholds
    Both,
}

/// Helper function to get the outlier detection method
///
/// ### Params
///
/// * `s` - String, type of test to run.
///
/// ### Returns
///
/// Option of the `OutlierDirection`
pub fn get_outlier_type(s: &str) -> Option<OutlierDirection> {
    match s.to_lowercase().as_str() {
        "below" => Some(OutlierDirection::Below),
        "above" => Some(OutlierDirection::Above),
        "twosided" => Some(OutlierDirection::Both),
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

    if len % 2 == 0 {
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

/// MAD outlier detection
///
/// ### Params
///
/// * `x` - Slice of values to check for outliers.
/// * `threshold` - Number of MADs to accept as not being an outlier.
/// * `direction` - Direction to check for outliers (below, above, or both).
///
/// ### Returns
///
/// A tuple with a vector of booleans indicating whether each value in the input
/// is an outlier and the value of applied margin.
pub fn mad_outlier(x: &[f64], threshold: f64, direction: OutlierDirection) -> (Vec<bool>, f64) {
    let median_val = match median(x) {
        Some(m) => m,
        None => return (vec![], 0_f64),
    };

    let mad_val = match mad(x) {
        Some(m) => m,
        None => return (vec![], 0_f64),
    };

    let margin = threshold * mad_val;

    let res = x
        .iter()
        .map(|&v| match direction {
            OutlierDirection::Below => v < median_val - margin,
            OutlierDirection::Above => v > median_val + margin,
            OutlierDirection::Both => (v - median_val).abs() > margin,
        })
        .collect::<Vec<bool>>();

    (res, margin)
}

/// Implementation of the trigamma function (second derivative of ln(gamma(x)))
///
/// ### Params
///
/// * `x` - The value for which to calculate the trigamma function.
///
/// ### Returns
///
/// The trigamma value for the given input.
pub fn trigamma(x: f64) -> f64 {
    let mut x = x;
    let mut result = 0.0;

    // For small x, use recurrence to get to a larger value
    if x <= 5.0 {
        // Each step through this loop applies the recurrence relation
        // trigamma(x) = trigamma(x+1) + 1/x^2
        while x < 5.0 {
            result += 1.0 / (x * x);
            x += 1.0;
        }
    }

    // For large x, use the asymptotic series
    // Main term: 1/x + 1/(2*x^2) + ...
    let xx = x * x;
    result += 1.0 / x + 1.0 / (2.0 * xx) + 1.0 / (6.0 * xx * x);

    // Add Bernoulli number terms for better precision
    let xxx = xx * x;
    result += -1.0 / (30.0 * xxx * x) + 1.0 / (42.0 * xxx * xx * x) - 1.0 / (30.0 * xxx * xxx * x);

    result
}

/// Calculate the critical value using bootstrap resampling
///
/// ### Params
///
/// * `values` - Slice of values to resample from.
/// * `sample_size` - Number of samples to draw in the bootstrap sample.
/// * `alpha` - The significance level for the critical value.
/// * `seed` - Random seed for reproducibility.
///
/// ### Returns
///
/// The critical value at the specified alpha level.
pub fn calculate_critval(values: &[f64], sample_size: usize, alpha: &f64, seed: usize) -> f64 {
    let mut rng = StdRng::seed_from_u64(seed as u64);
    let mut random_sample: Vec<f64> = (0..sample_size)
        .map(|_| {
            let index = rng.random_range(0..values.len());
            values[index]
        })
        .collect();
    random_sample.sort_by(|a, b| b.partial_cmp(a).unwrap());
    let index = (alpha * random_sample.len() as f64).ceil() as usize;
    random_sample[index + 1]
}

//////////////////
// Effect sizes //
//////////////////

/// Calculate the Hedge's g effect size and its standard error
///
/// ### Params
///
/// * `mean_a` - The mean values of group a.
/// * `mean_b` - The mean values of group b.
/// * `std_a` - The standard deviations of group a.
/// * `std_b` - The standard deviations of group b.
/// * `n_a` - Number of samples in a.
/// * `n_b` - Number of samples in b.
/// * `small_sample_correction` - Apply a small sample correction? Recommended
///   when `n_a` + `n_b` ≤ 35.
///
/// ### Returns
///
/// A tuple with the effect sizes being the first element, and the standard
/// errors the second element.
pub fn hedge_g_effect(
    mean_a: &[f64],
    mean_b: &[f64],
    std_a: &[f64],
    std_b: &[f64],
    n_a: usize,
    n_b: usize,
    small_sample_correction: bool,
) -> EffectSizeRes {
    assert_same_len!(mean_a, mean_b, std_a, std_b);

    let total_n = (n_a + n_b) as f64;
    let res: Vec<(f64, f64)> = mean_a
        .par_iter()
        .zip(mean_b.par_iter())
        .zip(std_a.par_iter())
        .zip(std_b.par_iter())
        .map(|(((mean_a, mean_b), std_a), std_b)| {
            let pooled_sd = (((n_a - 1) as f64 * std_a.powi(2) + (n_b - 1) as f64 * std_b.powi(2))
                / ((n_a + n_b - 2) as f64))
                .sqrt();
            let effect_size = (mean_a - mean_b) / pooled_sd;
            // Small sample correction if needed
            let effect_size = if small_sample_correction {
                let correction_factor =
                    ((total_n - 3.0) / (total_n - 2.25)) * ((total_n - 2.0) / total_n).sqrt();
                correction_factor * effect_size
            } else {
                effect_size
            };
            let standard_error = ((total_n / (n_a as f64 * n_b as f64))
                + (effect_size.powi(2) / (2.0 * total_n)))
                .sqrt();
            (effect_size, standard_error)
        })
        .collect();

    let mut effect_sizes: Vec<f64> = Vec::with_capacity(res.len());
    let mut standard_errors: Vec<f64> = Vec::with_capacity(res.len());

    for (effect_size, standard_error) in res {
        effect_sizes.push(effect_size);
        standard_errors.push(standard_error);
    }

    (effect_sizes, standard_errors)
}

////////////
// MANOVA //
////////////

/// ManovaResults
///
/// ### Fields
///
/// * `sscp_between` - Between-groups SSCP matrix
/// * `sscp_within` - Within-groups SSCP matrix
/// * `sscp_total` - Total SSCP matrix
/// * `df_between` - Degrees of freedom between groups
/// * `df_within` - Degrees of freedom within groups
/// * `df_total` - Total degrees of freedom
/// * `n_vars` - Number of variables
/// * `group_means` - Means for each group
/// * `overall_mean` - Overall means
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct ManovaResult {
    pub sscp_between: Mat<f64>,
    pub sscp_within: Mat<f64>,
    pub sscp_total: Mat<f64>,
    pub df_between: usize,
    pub df_within: usize,
    pub df_total: usize,
    pub n_vars: usize,
    pub group_means: Vec<Vec<f64>>,
    pub overall_mean: Vec<f64>,
}

#[allow(dead_code)]
impl ManovaResult {
    /// Calculate Wilks' Lambda statistic
    ///
    /// ### Returns
    ///
    /// Function will return the Wilks' Lambda
    pub fn wilks_lambda(&self) -> f64 {
        let w_plus_b = &self.sscp_within + &self.sscp_between;
        let det_w = self.sscp_within.determinant();
        let det_total = w_plus_b.determinant();
        det_w / det_total
    }

    /// Calculate Pillai's trace statistic
    ///
    /// ### Returns
    ///
    /// Function will return Pillai's trace
    pub fn pillai_trace(&self) -> f64 {
        let w_plus_b = &self.sscp_within + &self.sscp_between;
        let w_plus_b_pu = w_plus_b.partial_piv_lu();
        let w_plus_in = w_plus_b_pu.inverse();
        let h_times_inv = &self.sscp_between * w_plus_in;

        h_times_inv.diagonal().column_vector().iter().sum::<f64>()
    }

    /// Get F-statistic and p-value for Wilks' Lambda
    ///
    /// ### Returns
    ///
    /// A tuple with the F statistic and p-value according to Wilks'
    pub fn wilks_f_test(&self) -> (f64, f64) {
        let lambda = self.wilks_lambda();
        let p = self.n_vars as f64;
        let q = self.df_between as f64;
        let n = (self.df_within + self.df_between + 1) as f64;

        let t = ((p * p + q * q - 5.0).max(0.0)).sqrt();
        let w = n - (p + q + 1.0) / 2.0;
        let df1 = p * q;
        let df2 = w * t - (p * q - 2.0) / 2.0;

        let lambda_root = if t > 1.0 {
            lambda.powf(1.0 / t)
        } else {
            lambda
        };

        let f_stat = ((1.0 - lambda_root) / lambda_root) * (df2 / df1);

        let f_dist = FisherSnedecor::new(df1, df2).unwrap();
        let p_value = 1.0 - f_dist.cdf(f_stat);

        (f_stat, p_value)
    }

    /// Get F-statistic and p-value for Wilks' Lambda
    ///
    /// (Version used in R)
    ///
    /// ### Returns
    ///
    /// A tuple with the F statistic and p-value according to Pillai
    pub fn pillai_f_test(&self) -> (f64, f64) {
        let pillai = self.pillai_trace();
        let p = self.n_vars as f64;
        let q = self.df_between as f64;
        let n = self.df_within as f64;

        // F approximation for Pillai's trace
        let df1 = p * q;
        let df2 = q * (n - p + 1.0);

        let f_stat = (pillai / (q - pillai)) * (df2 / df1);

        let f_dist = FisherSnedecor::new(df1, df2).unwrap();
        let p_value = 1.0 - f_dist.cdf(f_stat);

        (f_stat, p_value)
    }
}

/// ManovaSummary
///
/// ### Fields
///
/// * `wilks_lambda` - Wilks' lambda value
/// * `pillai_trace` - Pillai's trace value
/// * `df_between` - Degrees of freedem between groups
/// * `df_within` - Degrees of freedom within groups
/// * `f_stat_wilk` - F statistic according to Wilk
/// * `p_val_wilk` - P-value according to Wilk
/// * `f_stat_pillai` - F statistic according to Pillai
/// * `p_val_pillai` - P-value according to Pillai
#[derive(Debug)]
#[allow(dead_code)]
pub struct ManovaSummary {
    pub wilks_lambda: f64,
    pub pillai_trace: f64,
    pub df_between: usize,
    pub df_within: usize,
    pub f_stat_wilk: f64,
    pub p_val_wilk: f64,
    pub f_stat_pillai: f64,
    pub p_val_pillai: f64,
}

impl ManovaSummary {
    /// Get the summary results from a `ManovaRes`
    ///
    /// ### Params
    ///
    /// * `res` - The calculated ManovaResults
    ///
    /// ### Returns
    ///
    /// The `ManovaSummary`.
    pub fn from_manova_res(res: &ManovaResult) -> Self {
        let (f_stat_wilk, p_val_wilk) = res.wilks_f_test();
        let (f_stat_pillai, p_val_pillai) = res.pillai_f_test();

        ManovaSummary {
            wilks_lambda: res.wilks_lambda(),
            pillai_trace: res.pillai_trace(),
            df_between: res.df_between,
            df_within: res.df_within,
            f_stat_wilk,
            p_val_wilk,
            f_stat_pillai,
            p_val_pillai,
        }
    }
}

///////////
// ANOVA //
///////////

/// AnovaSummary (based on MANOVA models)
///
/// ### Fields
///
/// * `variable_index` - Wilks' lambda value
/// * `ss_between` - Pillai's trace value
/// * `ss_within` - Degrees of freedem between groups
/// * `ms_between` - Degrees of freedom within groups
/// * `ms_within` - F statistic according to Wilk
/// * `f_stat` - P-value according to Wilk
/// * `p_val` - F statistic according to Pillai
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct AnovaSummary {
    pub variable_index: usize,
    pub ss_between: f64,
    pub ss_within: f64,
    pub ms_between: f64,
    pub ms_within: f64,
    pub f_stat: f64,
    pub p_val: f64,
}

pub fn summary_aov(res: &ManovaResult) -> Vec<AnovaSummary> {
    let mut aov_res = Vec::with_capacity(res.n_vars);

    for var_idx in 0..res.n_vars {
        let ss_between = res.sscp_between[(var_idx, var_idx)];
        let ss_within = res.sscp_within[(var_idx, var_idx)];
        let ms_between = ss_between / res.df_between as f64;
        let ms_within = ss_within / res.df_within as f64;
        let f_stat = ms_between / ms_within;

        let f_dist = FisherSnedecor::new(res.df_between as f64, res.df_within as f64).unwrap();
        let pval = 1.0 - f_dist.cdf(f_stat);

        aov_res.push(AnovaSummary {
            variable_index: var_idx,
            ss_between,
            ss_within,
            ms_between,
            ms_within,
            f_stat,
            p_val: pval,
        });
    }

    aov_res
}
