//! Hotfix for mitch at large scale

use bixverse_rs::core::math::matrix_helpers::*;
use bixverse_rs::core::math::vector_helpers::*;
use bixverse_rs::prelude::*;
use faer::{linalg::solvers::DenseSolveCore, Mat, MatRef};
use rayon::prelude::*;
use statrs::distribution::ContinuousCDF;
use statrs::distribution::FisherSnedecor;

// -- Statistical tests - key parts from stats.rs in bixverse-rs --

////////////
// MANOVA //
////////////

/// ManovaResults
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct ManovaResult<T>
where
    T: BixverseFloat,
{
    /// Between-groups SSCP matrix
    pub sscp_between: Mat<T>,
    /// Within-groups SSCP matrix
    pub sscp_within: Mat<T>,
    /// Total SSCP matrix
    pub sscp_total: Mat<T>,
    /// Degrees of freedom between groups
    pub df_between: usize,
    /// Degrees of freedom within groups
    pub df_within: usize,
    /// Total degrees of freedom
    pub df_total: usize,
    /// Number of variables
    pub n_vars: usize,
    /// Means for each group
    pub group_means: Vec<Vec<T>>,
    /// Overall means
    pub overall_mean: Vec<T>,
}

#[allow(dead_code)]
impl<T> ManovaResult<T>
where
    T: BixverseFloat + std::iter::Sum,
{
    /// Non-zero eigenvalue of E^-1 H for two-group MANOVA.
    ///
    /// For rank-1 H, the unique non-zero generalised eigenvalue equals
    /// trace(E^-1 H). All MANOVA statistics derive from this without
    /// computing det(E) / det(E+H), which overflows f64 for moderately
    /// large p (entries of E scale ~ n^3/12 for ranked data).
    fn rank1_eigenvalue(&self) -> T {
        debug_assert_eq!(
            self.df_between, 1,
            "rank1_eigenvalue assumes two-group MANOVA (df_between == 1)"
        );
        let e_inv = self.sscp_within.partial_piv_lu().inverse();
        let prod = &e_inv * &self.sscp_between;
        let trace: T = prod.diagonal().column_vector().iter().copied().sum::<T>();
        // Numerical noise can make this slightly negative for ill-conditioned E.
        trace.max(T::zero())
    }

    /// Wilks' Lambda for two-group MANOVA: Λ = 1 / (1 + λ).
    pub fn wilks_lambda(&self) -> T {
        let lambda = self.rank1_eigenvalue();
        T::one() / (T::one() + lambda)
    }

    /// Pillai's trace for two-group MANOVA: V = λ / (1 + λ).
    pub fn pillai_trace(&self) -> T {
        let lambda = self.rank1_eigenvalue();
        lambda / (T::one() + lambda)
    }

    /// Exact F-test for two-group MANOVA (Hotelling's T^2 form).
    ///
    /// F = ((n - p - 1) / p) * λ,   df = (p, n - p - 1)
    ///
    /// Returns (NaN, 1.0) when the inputs are degenerate
    /// (n - p - 1 <= 0, non-finite eigenvalue, etc.).
    fn two_group_f_test(&self) -> (T, T) {
        let lambda = self.rank1_eigenvalue();
        let p = T::from_usize(self.n_vars).unwrap();
        let n = T::from_usize(self.df_within + self.df_between + 1).unwrap();
        let df2 = n - p - T::one();

        if !lambda.is_finite() || df2 <= T::zero() {
            return (T::nan(), T::one());
        }

        let f_stat = (df2 / p) * lambda;

        let df1_f64 = p.to_f64().unwrap();
        let df2_f64 = df2.to_f64().unwrap();
        let f_f64 = f_stat.to_f64().unwrap();

        if !f_f64.is_finite() || f_f64 < 0.0 {
            return (f_stat, T::one());
        }

        let f_dist = match FisherSnedecor::new(df1_f64, df2_f64) {
            Ok(d) => d,
            Err(_) => return (f_stat, T::one()),
        };
        let p_value = T::from_f64(1.0 - f_dist.cdf(f_f64)).unwrap();
        (f_stat, p_value)
    }

    /// F-statistic and p-value derived from Wilks' Lambda.
    /// For two-group MANOVA this coincides with the Pillai F-test.
    pub fn wilks_f_test(&self) -> (T, T) {
        self.two_group_f_test()
    }

    /// F-statistic and p-value derived from Pillai's trace.
    /// For two-group MANOVA this coincides with the Wilks F-test.
    pub fn pillai_f_test(&self) -> (T, T) {
        self.two_group_f_test()
    }
}

/// ManovaSummary
#[derive(Debug)]
#[allow(dead_code)]
pub struct ManovaSummary<T>
where
    T: BixverseFloat,
{
    /// Wilks' lambda value
    pub wilks_lambda: T,
    /// Pillai's trace value
    pub pillai_trace: T,
    /// Degrees of freedom between groups
    pub df_between: usize,
    /// Degrees of freedom within groups
    pub df_within: usize,
    /// F statistic according to Wilk
    pub f_stat_wilk: T,
    /// P-value according to Wilk
    pub p_val_wilk: T,
    /// F statistic according to Pillai
    pub f_stat_pillai: T,
    /// P-value according to Pillai
    pub p_val_pillai: T,
}

impl<T> ManovaSummary<T>
where
    T: BixverseFloat + std::iter::Sum,
{
    /// Get the summary results from a `ManovaRes`
    ///
    /// ### Params
    ///
    /// * `res` - The calculated ManovaResults
    ///
    /// ### Returns
    ///
    /// The `ManovaSummary`.
    pub fn from_manova_res(res: &ManovaResult<T>) -> Self {
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
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct AnovaSummary<T>
where
    T: BixverseFloat,
{
    /// Variable index
    pub variable_index: usize,
    /// Sum of squares between groups
    pub ss_between: T,
    /// Sum of squares within groups
    pub ss_within: T,
    /// Mean square between groups
    pub ms_between: T,
    /// Mean square within groups
    pub ms_within: T,
    /// F statistic
    pub f_stat: T,
    /// P-value
    pub p_val: T,
}

/// Generates from MANOVE results the AnovaSummary
///
/// ### Params
///
/// * `res` - The MANOVA result to analyse
///
/// ### Returns
///
/// A vector of AnovaSummaries
pub fn summary_aov<T>(res: &ManovaResult<T>) -> Vec<AnovaSummary<T>>
where
    T: BixverseFloat,
{
    let mut aov_res = Vec::with_capacity(res.n_vars);
    let df_between_t = T::from_usize(res.df_between).unwrap();
    let df_within_t = T::from_usize(res.df_within).unwrap();

    let f_dist = FisherSnedecor::new(
        df_between_t.to_f64().unwrap(),
        df_within_t.to_f64().unwrap(),
    )
    .ok();

    for var_idx in 0..res.n_vars {
        let ss_between = res.sscp_between[(var_idx, var_idx)];
        let ss_within = res.sscp_within[(var_idx, var_idx)];
        let ms_between = ss_between / df_between_t;
        let ms_within = ss_within / df_within_t;

        let (f_stat, p_val) = if ms_within <= T::zero() {
            (T::nan(), T::one())
        } else {
            let f = ms_between / ms_within;
            let f_f64 = f.to_f64().unwrap();
            let pv = match (&f_dist, f_f64.is_finite() && f_f64 >= 0.0) {
                (Some(d), true) => T::from_f64(1.0 - d.cdf(f_f64)).unwrap(),
                _ => T::one(),
            };
            (f, pv)
        };

        aov_res.push(AnovaSummary {
            variable_index: var_idx,
            ss_between,
            ss_within,
            ms_between,
            ms_within,
            f_stat,
            p_val,
        });
    }

    aov_res
}

// -- Statistical tests - key parts from mitch.rs in bixverse-rs --

//////////////////
// Type aliases //
//////////////////

/// Type alias for processed MitchPathways
///
/// ### Slots
///
/// * `0` - Name of the pathway
/// * `1` - Index positions of the pathway relative to the input matrix
pub type MitchPathways = (Vec<String>, Vec<Vec<usize>>);

////////////////////
// Result structs //
////////////////////

/// Structure to store Mitch results
///
/// ### Fields
///
/// * `pathway_name` - Name of the pathway
/// * `pathway_size` - Size of the pathway
/// * `manova_pval` - P-value of the MANOVA test on the ranked data
/// * `scores` - The scores for each tested contrast for this pathway
/// * `anova_pvals` - The p-values for the individual contrasts based on the
///   ANOVA on top of the MANOVA results
/// * `s_dist` - Calculated distances from the hypotenuse.
/// * `mysd` - The standard deviation of the scores for this pathway.
#[derive(Clone, Debug)]
pub struct MitchResult<'a, T> {
    /// Name of the pathway
    pub pathway_name: &'a str,
    /// Size of the pathway
    pub pathway_size: usize,
    /// P-value of the MANOVA test on the ranked data
    pub manova_pval: T,
    /// The scores for each tested contrast for this pathway
    pub scores: Vec<T>,
    /// The p-values for the individual contrasts based on the ANOVA on top of
    /// the MANOVA results
    pub anova_pvals: Vec<T>,
    /// Calculated distances from the hypotenuse.
    pub s_dist: T,
    /// The standard deviation of the scores for this pathway.
    pub mysd: T,
}

///////////////
// Functions //
///////////////

/// Calculate the MANOVA results for a given pathway
///
/// ### Params
///
/// * `x` - The pre-ranked matrix.
/// * `group1_indices` - The row index positions for which genes belong to the
///   pathway.
///
/// ### Return
///
/// Returns the MANOVA results for this pathway.
pub fn manova_mitch<T: BixverseFloat>(x: MatRef<T>, group1_indices: &[usize]) -> ManovaResult<T> {
    let (n, p) = x.shape();

    assert!(
        group1_indices.iter().all(|&idx| idx < n),
        "All indices must be less than the number of rows"
    );

    let mut is_group1 = vec![false; n];
    for &idx in group1_indices {
        is_group1[idx] = true;
    }

    let group0_indices: Vec<usize> = (0..n).filter(|&i| !is_group1[i]).collect();

    let n1 = group1_indices.len();
    let n0 = group0_indices.len();

    let mut x0 = Mat::zeros(n0, p);
    let mut x1 = Mat::zeros(n1, p);

    for (new_i, &orig_i) in group0_indices.iter().enumerate() {
        for j in 0..p {
            x0[(new_i, j)] = x[(orig_i, j)];
        }
    }

    for (new_i, &orig_i) in group1_indices.iter().enumerate() {
        for j in 0..p {
            x1[(new_i, j)] = x[(orig_i, j)];
        }
    }

    let mean_0 = col_means(x0.as_ref());
    let mean_1 = col_means(x1.as_ref());
    let mean_overall = col_means(x);

    let x_centered = scale_matrix_col(&x.as_ref(), false);
    let x0_centered = scale_matrix_col(&x0.as_ref(), false);
    let x1_centered = scale_matrix_col(&x1.as_ref(), false);

    let sscp_within =
        x0_centered.transpose() * &x0_centered + x1_centered.transpose() * &x1_centered;
    let sscp_total = x_centered.transpose() * &x_centered;
    let sscp_between = &sscp_total - &sscp_within;

    ManovaResult {
        sscp_between,
        sscp_within,
        sscp_total,
        df_between: 1,
        df_within: n - 2,
        df_total: n - 1,
        n_vars: p,
        group_means: vec![mean_0, mean_1],
        overall_mean: mean_overall,
    }
}

/// Calculates the mitch-specific ranks for a given contrast based on the tied
/// method
///
/// ### Params
///
/// * `mat` - The matrix to rank
///
/// ### Returns
///
/// The mitched-ranked matrix
pub fn mitch_rank<T: BixverseFloat>(mat: &MatRef<T>) -> Mat<T> {
    let mut ranked_mat = Mat::zeros(mat.nrows(), mat.ncols());

    ranked_mat
        .par_col_iter_mut()
        .enumerate()
        .for_each(|(col_idx, mut col)| {
            let original_col: Vec<T> = mat.col(col_idx).iter().copied().collect();
            let mut zeros = T::zero();
            let mut neg = T::zero();
            for x in &original_col {
                if *x == T::zero() {
                    zeros += T::one();
                } else if *x < T::zero() {
                    neg += T::one();
                }
            }
            let adj = neg + (zeros / T::from_f64(2.0).unwrap());

            let ranks = rank_vector(&original_col);
            let ranks = ranks.iter().map(|x| *x - adj).collect::<Vec<T>>();

            for (row_idx, &rank) in ranks.iter().enumerate() {
                col[row_idx] = rank;
            }
        });

    ranked_mat
}

/// Wrapper function to process a given pathway
///
/// ### Params
///
/// * `ranked_mat` - The ranked matrix
/// * `pathway_name` - Name of the pathway/gene set that is being tested.
/// * `pathway_indices` - Index positions which genes (rows) belong to this given
///   pathway
///
/// ### Returns
///
/// A `MitchResult` structure with the results for this pathway.
pub fn process_mitch_pathway<'a, T>(
    ranked_mat: MatRef<T>,
    pathway_name: &'a str,
    pathway_indices: &[usize],
) -> MitchResult<'a, T>
where
    T: BixverseFloat + std::iter::Sum,
{
    let nrow = T::from_usize(ranked_mat.nrows()).unwrap();
    let manova_res: ManovaResult<T> = manova_mitch(ranked_mat, pathway_indices);
    let sum_manova: ManovaSummary<T> = ManovaSummary::from_manova_res(&manova_res);
    let sum_anova = summary_aov(&manova_res);

    let p_manova = sum_manova.p_val_pillai;
    let mut p_aovs = Vec::with_capacity(sum_anova.len());

    for sum in sum_anova {
        p_aovs.push(sum.p_val);
    }

    let scores = &manova_res.group_means[0]
        .iter()
        .zip(&manova_res.group_means[1])
        .map(|(mean_0, mean_1)| (T::from_f64(2.0).unwrap() * (*mean_1 - *mean_0)) / nrow)
        .collect::<Vec<T>>();

    let sd_scores = standard_deviation(scores);
    let hypotenuse = scores
        .iter()
        .map(|x| x.powi(2))
        .fold(T::zero(), |acc, x| acc + x)
        .sqrt();

    MitchResult {
        pathway_name,
        pathway_size: pathway_indices.len(),
        manova_pval: p_manova,
        scores: scores.clone(),
        anova_pvals: p_aovs,
        s_dist: hypotenuse,
        mysd: sd_scores,
    }
}

///////////
// Tests //
///////////

#[test]
fn test_manova_two_group_analytic() {
    // Group 0: rows 0, 1; Group 1: rows 2, 3
    // x = [[1, 4], [2, 1], [5, 7], [6, 3]]
    //
    // Hand-computed:
    //   sscp_within = [[1.0, -3.5], [-3.5, 12.5]]
    //   d = mean_1 - mean_0 = [4.0, 2.5]
    //   c = n1*n2/n = 1
    //   lambda = c * d' * E^-1 * d = 1105
    //   Wilks = 1/1106, Pillai = 1105/1106, Wilks + Pillai = 1
    //   F = ((n - p - 1)/p) * lambda = (1/2)*1105 = 552.5, df = (2, 1)
    let mat = Mat::from_fn(4, 2, |i, j| match (i, j) {
        (0, 0) => 1.0f64,
        (0, 1) => 4.0,
        (1, 0) => 2.0,
        (1, 1) => 1.0,
        (2, 0) => 5.0,
        (2, 1) => 7.0,
        (3, 0) => 6.0,
        (3, 1) => 3.0,
        _ => unreachable!(),
    });
    let res = manova_mitch(mat.as_ref(), &[2, 3]);

    let wilks = res.wilks_lambda();
    let pillai = res.pillai_trace();
    let recovered_lambda = (1.0 - wilks) / wilks;

    assert!((recovered_lambda - 1105.0).abs() < 1e-6);
    assert!((wilks + pillai - 1.0).abs() < 1e-12);
    assert!((wilks - 1.0 / 1106.0).abs() < 1e-9);
    assert!((pillai - 1105.0 / 1106.0).abs() < 1e-9);

    let (f_w, _) = res.wilks_f_test();
    let (f_p, p_p) = res.pillai_f_test();
    assert!((f_w - 552.5).abs() < 1e-6);
    assert!((f_p - 552.5).abs() < 1e-6);
    assert!((0.0..=1.0).contains(&p_p));
}

#[test]
fn test_manova_no_overflow_large_p() {
    // Regression test: previously panicked via det() overflow → NaN →
    // FisherSnedecor::cdf XOutOfRange. Mimics Mitch ranked-data magnitudes
    // (n ~ 3000, sscp diagonals ~ n^3/12 ~ 2e9).
    let p = 80;
    let n_total = 2949usize;
    let sigma = (n_total as f64).powi(3) / 12.0;

    let sscp_within = Mat::from_fn(p, p, |i, j| if i == j { sigma } else { sigma * 0.3 });
    // Small between-group signal; H entries ~ sigma * 1e-3
    let sscp_between = Mat::from_fn(
        p,
        p,
        |i, j| {
            if i == j {
                sigma * 1.5e-3
            } else {
                sigma * 5e-4
            }
        },
    );
    let sscp_total = &sscp_within + &sscp_between;

    let res = ManovaResult::<f64> {
        sscp_between,
        sscp_within,
        sscp_total,
        df_between: 1,
        df_within: n_total - 2,
        df_total: n_total - 1,
        n_vars: p,
        group_means: vec![vec![0.0; p], vec![0.0; p]],
        overall_mean: vec![0.0; p],
    };

    let wilks = res.wilks_lambda();
    let pillai = res.pillai_trace();
    let (f_w, p_w) = res.wilks_f_test();
    let (f_p, p_p) = res.pillai_f_test();

    assert!(wilks.is_finite() && wilks > 0.0 && wilks <= 1.0);
    assert!((0.0..1.0).contains(&pillai) && pillai.is_finite());
    assert!((wilks + pillai - 1.0).abs() < 1e-9);
    assert!(f_w.is_finite() && f_w >= 0.0);
    assert!(f_p.is_finite() && f_p >= 0.0);
    assert!((0.0..=1.0).contains(&p_w));
    assert!((0.0..=1.0).contains(&p_p));
    // Wilks and Pillai F-tests are identical for two groups.
    assert!((f_w - f_p).abs() < 1e-9);
    assert!((p_w - p_p).abs() < 1e-9);
}

#[test]
fn test_summary_aov_constant_column() {
    // ss_within == 0 used to feed NaN into FisherSnedecor and panic.
    let p = 2;
    let sscp_within = Mat::from_fn(p, p, |i, j| {
        if i == 0 && j == 0 {
            0.0
        } else if i == j {
            10.0
        } else {
            0.0
        }
    });
    let sscp_between = Mat::from_fn(p, p, |i, j| if i == j { 5.0 } else { 0.0 });
    let sscp_total = &sscp_within + &sscp_between;

    let res = ManovaResult::<f64> {
        sscp_between,
        sscp_within,
        sscp_total,
        df_between: 1,
        df_within: 10,
        df_total: 11,
        n_vars: p,
        group_means: vec![vec![0.0; p], vec![0.0; p]],
        overall_mean: vec![0.0; p],
    };

    let aov = summary_aov(&res);
    assert_eq!(aov.len(), 2);
    assert!(aov[0].f_stat.is_nan());
    assert!((aov[0].p_val - 1.0).abs() < 1e-12);
    assert!(aov[1].f_stat.is_finite() && aov[1].f_stat > 0.0);
    assert!((0.0..=1.0).contains(&aov[1].p_val));
}
