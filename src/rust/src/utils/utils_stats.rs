use faer::{ColRef, Mat, MatRef};
use rand::prelude::*;
use rand::seq::SliceRandom;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use statrs::distribution::{Continuous, ContinuousCDF, Normal};
use std::ops::{Add, Div};

use crate::assert_same_len;
use crate::helpers::linalg::col_sums;

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

///////////////////////
// Generic functions //
///////////////////////

/// Split a vector randomly into two chunks
///
/// Splits a vector randomly into two of [..x] and the other [x..]
///
/// ### Params
///
/// * `vec` - Slice of the vector you want to split
/// * `x` - Length of the first vector; the rest will be put into the second vector
/// * `seed` - Seed for reproducibility
///
/// ### Returns
///
/// A tuple of the pieces of the vector
pub fn split_vector_randomly(vec: &[f64], x: usize, seed: u64) -> (Vec<f64>, Vec<f64>) {
    let mut rng = StdRng::seed_from_u64(seed);
    let mut shuffled = vec.to_vec();
    shuffled.shuffle(&mut rng);

    let (first_set, second_set) = shuffled.split_at(x);

    (first_set.to_vec(), second_set.to_vec())
}

/// Calculate the set similarity.
///
/// ### Params
///
/// * `s_1` - The first HashSet.
/// * `s_2` - The second HashSet.
/// * `overlap_coefficient` - Shall the overlap coefficient be returned or the
///   Jaccard similarity
///
/// ### Return
///
/// The Jaccard similarity or overlap coefficient.
pub fn set_similarity(
    s_1: &FxHashSet<&String>,
    s_2: &FxHashSet<&String>,
    overlap_coefficient: bool,
) -> f64 {
    let i = s_1.intersection(s_2).count() as u64;
    let u = if overlap_coefficient {
        std::cmp::min(s_1.len(), s_2.len()) as u64
    } else {
        s_1.union(s_2).count() as u64
    };
    i as f64 / u as f64
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

//////////////////////
// Vector functions //
//////////////////////

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

////////////////////////////
// Radial Basis functions //
////////////////////////////

/// Enum for the RBF function
#[derive(Debug)]
pub enum RbfType {
    Gaussian,
    Bump,
    InverseQuadratic,
}

/// Parsing the RBF function
///
/// ### Params
///
/// * `s` - String to transform into `RbfType`
///
/// ### Returns
///
/// Returns the `RbfType`
pub fn parse_rbf_types(s: &str) -> Option<RbfType> {
    match s.to_lowercase().as_str() {
        "gaussian" => Some(RbfType::Gaussian),
        "bump" => Some(RbfType::Bump),
        "inverse_quadratic" => Some(RbfType::InverseQuadratic),
        _ => None,
    }
}

/// Gaussian Radial Basis function
///
/// Applies a Gaussian Radial Basis function on a vector of distances with the
/// following formula:
/// ```
/// φ(r) = e^(-(εr)²)
/// ```
///
/// ### Params
///
/// * `dist` - Vector of distances
/// * `epsilon` - Shape parameter controlling function width
///
/// ### Returns
///
/// The resulting affinity vector
pub fn rbf_gaussian(dist: &[f64], epsilon: &f64) -> Vec<f64> {
    dist.par_iter()
        .map(|x| f64::exp(-((x * *epsilon).powi(2))))
        .collect()
}

/// Gaussian Radial Basis function for matrices.
///
/// Applies a Gaussian Radial Basis function on a matrix of distances with the
/// following formula:
///
/// ```
/// φ(r) = e^(-(εr)²)
/// ```
///
/// ### Params
///
/// * `dist` - Matrix of distances
/// * `epsilon` - Shape parameter controlling function width
///
/// ### Returns
///
/// The affinity matrix
pub fn rbf_gaussian_mat(dist: MatRef<f64>, epsilon: &f64) -> Mat<f64> {
    let ncol = dist.ncols();
    let nrow = dist.nrows();
    Mat::from_fn(nrow, ncol, |i, j| {
        let x = dist.get(i, j);
        f64::exp(-((x * epsilon).powi(2)))
    })
}

/// Bump Radial Basis function
///
/// Applies a Bump Radial Basis function on a vector of distances with the
/// following formula:
/// ```
/// φ(r) = { exp(-1/(1-(εr)²)) + 1,  if εr < 1
///        { 0,                      if εr ≥ 1
/// ```
///
/// ### Params
///
/// * `dist` - Vector of distances
/// * `epsilon` - Shape parameter controlling function width
///
/// ### Returns
///
/// The resulting affinity vector
pub fn rbf_bump(dist: &[f64], epsilon: &f64) -> Vec<f64> {
    dist.par_iter()
        .map(|x| {
            if *x < (1.0 / epsilon) {
                f64::exp(-(1_f64 / (1_f64 - (*epsilon * x).powi(2))) + 1_f64)
            } else {
                0_f64
            }
        })
        .collect()
}

/// Bump Radial Basis function for matrices
///
/// Applies a Bump Radial Basis function on a matrix of distances with the
/// following formula:
/// ```
/// φ(r) = { exp(-1/(1-(εr)²)) + 1,  if εr < 1
///        { 0,                      if εr ≥ 1
/// ```
///
/// ### Params
///
/// * `dist` - Matrix of distances
/// * `epsilon` - Shape parameter controlling function width
///
/// ### Returns
///
/// The resulting affinity matrix
pub fn rbf_bump_mat(dist: MatRef<f64>, epsilon: &f64) -> Mat<f64> {
    let ncol = dist.ncols();
    let nrow = dist.nrows();
    Mat::from_fn(nrow, ncol, |i, j| {
        let x = dist.get(i, j);
        if *x < (1.0 / epsilon) {
            f64::exp(-(1_f64 / (1_f64 - (*epsilon * x).powi(2))) + 1_f64)
        } else {
            0_f64
        }
    })
}

/// Inverse quadratic RBF
///
/// Applies a Inverse Quadratic Radial Basis function on a vector of distances
/// with the following formula:
/// ```
/// φ(r) = 1/(1 + (εr)²)
/// ```
///
/// ### Params
///
/// * `dist` - Vector of distances
/// * `epsilon` - Shape parameter controlling function width
///
/// ### Return
///
/// The resulting affinity vector
pub fn rbf_inverse_quadratic(dist: &[f64], epsilon: &f64) -> Vec<f64> {
    dist.par_iter()
        .map(|x| 1.0 / (1.0 + (*epsilon * x).powi(2)))
        .collect()
}

/// Inverse quadratic RBF for matrices
///
/// Applies a Inverse Quadratic Radial Basis function on a matrix of distances
/// with the following formula:
/// ```
/// φ(r) = 1/(1 + (εr)²)
/// ```
///
/// ### Params
///
/// * `dist` - Matrix of distances
/// * `epsilon` - Shape parameter controlling function width
///
/// ### Returns
///
/// The resulting affinity matrix
pub fn rbf_inverse_quadratic_mat(dist: MatRef<f64>, epsilon: &f64) -> Mat<f64> {
    let ncol = dist.ncols();
    let nrow = dist.nrows();
    Mat::from_fn(nrow, ncol, |i, j| {
        let x = dist.get(i, j);
        1.0 / (1.0 + (*epsilon * x).powi(2))
    })
}

/////////////////////////////////
// Topological overlap measure //
/////////////////////////////////

/// Enum for the TOM function
#[derive(Debug)]
pub enum TomType {
    Version1,
    Version2,
}

/// Parsing the TOM type
pub fn parse_tom_types(s: &str) -> Option<TomType> {
    match s.to_lowercase().as_str() {
        "v1" => Some(TomType::Version1),
        "v2" => Some(TomType::Version2),
        _ => None,
    }
}

/// Calculates the topological overlap measure (TOM) for a given affinity matrix
///
/// The TOM quantifies the relative interconnectedness of two nodes in a network by measuring
/// how much they share neighbors relative to their connectivity. Higher TOM values indicate
/// nodes that are part of the same module or cluster.
///
/// ### Params
///
/// * `affinity_mat` - Symmetric affinity/adjacency matrix
/// * `signed` - Whether to use signed (absolute values) or unsigned connectivity
/// * `tom_type` - Algorithm version (Version1 or Version2)
///
/// ### Returns
/// Symmetric TOM matrix with values in [0,1] representing topological overlap
///
/// ### Mathematical Formulation
///
/// #### Connectivity
/// For node i: k_i = Σ_j |a_ij| (signed) or k_i = Σ_j a_ij (unsigned)
///
/// #### Shared Neighbors
/// For nodes i,j: l_ij = Σ_k (a_ik * a_kj) - a_ii*a_ij - a_ij*a_jj
///
/// #### TOM Calculation
///
/// **Version 1:**
/// - Numerator: a_ij + l_ij
/// - Denominator (unsigned): min(k_i, k_j) + 1 - a_ij
/// - Denominator (signed): min(k_i, k_j) + 1 - a_ij (if a_ij ≥ 0) or min(k_i, k_j) + 1 + a_ij (if a_ij < 0)
/// - TOM_ij = numerator / denominator
///
/// **Version 2:**
/// - Divisor (unsigned): min(k_i, k_j) + a_ij
/// - Divisor (signed): min(k_i, k_j) + a_ij (if a_ij ≥ 0) or min(k_i, k_j) - a_ij (if a_ij < 0)
/// - TOM_ij = 0.5 * (a_ij + l_ij/divisor)
pub fn calc_tom(affinity_mat: MatRef<f64>, signed: bool, tom_type: TomType) -> Mat<f64> {
    let n = affinity_mat.nrows();
    let mut tom_mat = Mat::<f64>::zeros(n, n);

    let connectivity = if signed {
        (0..n)
            .map(|i| (0..n).map(|j| affinity_mat.get(i, j).abs()).sum())
            .collect::<Vec<f64>>()
    } else {
        col_sums(affinity_mat.as_ref())
    };

    // Pre-compute for speed-ups -> Massive difference and good call @Claude
    let dot_products = affinity_mat.as_ref() * affinity_mat.as_ref();

    for i in 0..n {
        // Only upper triangle for speed...
        for j in (i + 1)..n {
            let a_ij = affinity_mat.get(i, j);

            // shared_neighbors = dot_product - excluded terms (k=i and k=j)
            let shared_neighbours = dot_products.get(i, j)
                - affinity_mat.get(i, i) * affinity_mat.get(i, j)
                - affinity_mat.get(i, j) * affinity_mat.get(j, j);

            let f_ki_kj = connectivity[i].min(connectivity[j]);

            let tom_value = match tom_type {
                TomType::Version1 => {
                    let numerator = a_ij + shared_neighbours;
                    let denominator = if signed {
                        if *a_ij >= 0.0 {
                            f_ki_kj + 1.0 - a_ij
                        } else {
                            f_ki_kj + 1.0 + a_ij
                        }
                    } else {
                        f_ki_kj + 1.0 - a_ij
                    };
                    numerator / denominator
                }
                TomType::Version2 => {
                    let divisor = if signed {
                        if *a_ij >= 0.0 {
                            f_ki_kj + a_ij
                        } else {
                            f_ki_kj - a_ij
                        }
                    } else {
                        f_ki_kj + a_ij
                    };
                    let neighbours = shared_neighbours / divisor;
                    0.5 * (a_ij + neighbours)
                }
            };

            tom_mat[(i, j)] = tom_value;
            tom_mat[(j, i)] = tom_value;
        }
    }

    tom_mat
}

///////////////////////
// Cluster stability //
///////////////////////

/// Helper functions to calculate the intersection of sorted usize vectors
///
/// ### Params
///
/// * `a` - (Sorted) slice of usize
/// * `b` - (Sorted) slice of usize
///
/// ### Returns
///
/// Intersection between the two sorted slices.
fn intersection_size_sorted(a: &[usize], b: &[usize]) -> usize {
    let mut count = 0;
    let mut i = 0;
    let mut j = 0;

    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Equal => {
                count += 1;
                i += 1;
                j += 1;
            }
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
        }
    }

    count
}

/// Function that assesses the cluster stability.
///
/// ### Params
///
/// * `cluster_matrix` - A matrix with the columns representing the bootstraps,
///   the rows the features and the values which cluster the feature belongs to.
///
/// ### Returns
///
/// Returns tuple of (average Jaccard similarities, standard deviations of the
/// Jaccard similarities).
pub fn cluster_stability(cluster_matrix: &MatRef<i32>) -> Vec<(f64, f64)> {
    let n_features = cluster_matrix.nrows();
    let n_iter = cluster_matrix.ncols();

    // Pre-compute cluster membership maps for all bootstraps
    let bootstrap_cluster_maps: Vec<FxHashMap<i32, Vec<usize>>> = (0..n_iter)
        .into_par_iter()
        .map(|boot_idx| {
            let mut clusters_map: FxHashMap<i32, Vec<usize>> = FxHashMap::default();
            for feature_idx in 0..n_features {
                let cluster_id = cluster_matrix[(feature_idx, boot_idx)];
                clusters_map
                    .entry(cluster_id)
                    .or_default()
                    .push(feature_idx);
            }
            clusters_map
        })
        .collect();

    // Process features in parallel with optimized memory usage
    (0..n_features)
        .into_par_iter()
        .map(|feature_idx| {
            let n_pairs = (n_iter * (n_iter - 1)) / 2;
            let mut jaccard_scores = Vec::with_capacity(n_pairs);

            for i in 0..(n_iter - 1) {
                for j in (i + 1)..n_iter {
                    let cluster_i = cluster_matrix[(feature_idx, i)];
                    let cluster_j = cluster_matrix[(feature_idx, j)];

                    let members_i = bootstrap_cluster_maps[i]
                        .get(&cluster_i)
                        .map(|v| v.as_slice())
                        .unwrap_or(&[]);

                    let members_j = bootstrap_cluster_maps[j]
                        .get(&cluster_j)
                        .map(|v| v.as_slice())
                        .unwrap_or(&[]);

                    let intersection_size = intersection_size_sorted(members_i, members_j);
                    let union_size = members_i.len() + members_j.len() - intersection_size;

                    let jaccard = if union_size == 0 {
                        0.0
                    } else {
                        intersection_size as f64 / union_size as f64
                    };
                    jaccard_scores.push(jaccard);
                }
            }

            let mean_jaccard = jaccard_scores.iter().sum::<f64>() / jaccard_scores.len() as f64;
            let variance = jaccard_scores
                .iter()
                .map(|x| (x - mean_jaccard).powi(2))
                .sum::<f64>()
                / jaccard_scores.len() as f64;
            let std_jaccard = variance.sqrt();

            (mean_jaccard, std_jaccard)
        })
        .collect()
}

/////////////////////////////////////
// Additional stat implementations //
/////////////////////////////////////

/// Implementation of the trigamma function (second derivative of ln(gamma(x)))
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

/// Transform Z-scores into p-values (assuming normality).
///
/// ### Params
///
/// * `z_scores` - The Z scores to transform to p-values
///
/// ### Returns
///
/// The p-value vector based on the Z scores (two sided)
pub fn z_scores_to_pval(z_scores: &[f64]) -> Vec<f64> {
    let normal = Normal::new(0.0, 1.0).unwrap();
    z_scores
        .iter()
        .map(|&z| {
            let abs_z = z.abs();
            if abs_z > 6.0 {
                // Deal with numeric precision problems for very large z-scores.
                let pdf = normal.pdf(abs_z);
                let p = pdf / abs_z * (1.0 - 1.0 / (abs_z * abs_z));
                2.0 * p
            } else {
                2.0 * (1.0 - normal.cdf(abs_z))
            }
        })
        .collect()
}

/// Calculates the critical value
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

////////////////////////
// Information theory //
////////////////////////

/// Calculate the mutual information between two column references
///
/// The columns need to be binned
///
/// ### Params
///
/// * `col_i` - Column reference to the first column to compare
/// * `col_j` - Column reference to the second column to compare
/// * `n_bins` - Optional number of bins. If not provided, will default to
///   `sqrt(nrows)`.
///
/// ### Returns
///
/// The mutual information between the two columns
pub fn calculate_mi(col_i: ColRef<usize>, col_j: ColRef<usize>, n_bins: Option<usize>) -> f64 {
    let n_rows = col_i.nrows();
    let n_bins = n_bins.unwrap_or_else(|| (n_rows as f64).sqrt() as usize);

    let mut joint_counts = vec![vec![0usize; n_bins]; n_bins];
    let mut marginal_i = vec![0usize; n_bins];
    let mut marginal_j = vec![0usize; n_bins];

    for i in 0..n_rows {
        let bin_i = col_i[i];
        let bin_j = col_j[i];

        joint_counts[bin_i][bin_j] += 1;
        marginal_i[bin_i] += 1;
        marginal_j[bin_j] += 1;
    }

    let n = n_rows as f64;
    let mut mi = 0.0;

    #[allow(clippy::needless_range_loop)]
    for i in 0..n_bins {
        for j in 0..n_bins {
            let joint_prob = joint_counts[i][j] as f64 / n;
            if joint_prob > 0.0 {
                let marginal_prob_i = marginal_i[i] as f64 / n;
                let marginal_prob_j = marginal_j[j] as f64 / n;

                mi += joint_prob * (joint_prob / (marginal_prob_i * marginal_prob_j)).ln();
            }
        }
    }

    mi
}

/// Calculates the joint entropy between two column references
///
/// ### Params
///
/// * `col_i` - Column reference to the first column to compare
/// * `col_j` - Column reference to the second column to compare
/// * `n_bins` - Optional number of bins. If not provided, will default to
///   `sqrt(nrows)`.
///
/// ### Returns
///
/// The joint entropy
pub fn calculate_joint_entropy(
    col_i: ColRef<usize>,
    col_j: ColRef<usize>,
    n_bins: Option<usize>,
) -> f64 {
    let n_rows = col_i.nrows();
    let n_bins = n_bins.unwrap_or_else(|| (n_rows as f64).sqrt() as usize);

    let mut joint_counts = vec![vec![0usize; n_bins]; n_bins];

    for i in 0..n_rows {
        joint_counts[col_i[i]][col_j[i]] += 1;
    }

    let n = n_rows as f64;
    let mut joint_entropy = 0.0;

    #[allow(clippy::needless_range_loop)]
    for i in 0..n_bins {
        for j in 0..n_bins {
            let joint_prob = joint_counts[i][j] as f64 / n;
            if joint_prob > 0.0 {
                joint_entropy -= joint_prob * joint_prob.ln();
            }
        }
    }

    joint_entropy
}

/// Calculates the entropy of a column reference
///
/// The column needs to be binned
///
/// ### Params
///
/// * `col` - Column reference for which to calculate entropy
/// * `n_bins` - Optional number of bins. If not provided, will default to
///   `sqrt(nrows)`.
///
/// ### Returns
///
/// The entropy of the column
pub fn calculate_entropy(col: ColRef<usize>, n_bins: Option<usize>) -> f64 {
    let n_rows = col.nrows();
    let n_bins = n_bins.unwrap_or_else(|| (n_rows as f64).sqrt() as usize);

    let mut counts = vec![0usize; n_bins];

    for i in 0..n_rows {
        counts[col[i]] += 1;
    }

    let n = n_rows as f64;
    let mut entropy = 0.0;

    for &count in &counts {
        if count > 0 {
            let prob = count as f64 / n;
            entropy -= prob * prob.ln();
        }
    }

    entropy
}
