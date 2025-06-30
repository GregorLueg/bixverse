use faer::{Mat, MatRef};
use rand::prelude::*;
use rand::seq::SliceRandom;
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use statrs::distribution::{Continuous, ContinuousCDF, Normal};

use crate::helpers_linalg::col_sums;

///////////
// Types //
///////////

/// A type alias that can be returned by the par_iter() functions.
pub type EffectSizeRes = (Vec<f64>, Vec<f64>);

///////////////////////
// Generic functions //
///////////////////////

/// Split a vector randomly into two chunks with one being [..x] and the other [x..]
pub fn split_vector_randomly(vec: &[f64], x: usize, seed: u64) -> (Vec<f64>, Vec<f64>) {
    let mut rng = StdRng::seed_from_u64(seed);
    let mut shuffled = vec.to_vec();
    shuffled.shuffle(&mut rng);

    let (first_set, second_set) = shuffled.split_at(x);

    (first_set.to_vec(), second_set.to_vec())
}

/// Calculate the set similarity. Options are Jaccard (similarity_index = False)
/// or the similarity index calculation.
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

//////////////////////
// Vector functions //
//////////////////////

/// Get the median of a vector
pub fn median(x: &[f64]) -> Option<f64> {
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
        Some((*median1 + *median2) / 2.0)
    } else {
        let (_, median, _) = data.select_nth_unstable_by(len / 2, |a, b| a.partial_cmp(b).unwrap());
        Some(*median)
    }
}

/// Calculate the median absolute deviation of a Vector
pub fn mad(data: &[f64]) -> Option<f64> {
    if data.is_empty() {
        return None;
    }

    let median_val = median(data)?; // Early return if median is None
    let deviations: Vec<f64> = data.iter().map(|&x| (x - median_val).abs()).collect();
    median(&deviations)
}

//////////////////
// Effect sizes //
//////////////////

/// Calculate the Hedge's g effect size and its standard error
pub fn hedge_g_effect(
    mean_a: &[f64],
    mean_b: &[f64],
    std_a: &[f64],
    std_b: &[f64],
    n_a: usize,
    n_b: usize,
    small_sample_correction: bool,
) -> EffectSizeRes {
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
pub fn parse_rbf_types(s: &str) -> Option<RbfType> {
    match s.to_lowercase().as_str() {
        "gaussian" => Some(RbfType::Gaussian),
        "bump" => Some(RbfType::Bump),
        "inverse_quadratic" => Some(RbfType::InverseQuadratic),
        _ => None,
    }
}

/// Gaussian Radial Basis function for vectors
pub fn rbf_gaussian(dist: &[f64], epsilon: &f64) -> Vec<f64> {
    dist.par_iter()
        .map(|x| f64::exp(-((x * *epsilon).powi(2))))
        .collect()
}

/// Gaussian Radial Basis function for matrices.
pub fn rbf_gaussian_mat(dist: MatRef<'_, f64>, epsilon: &f64) -> Mat<f64> {
    let ncol = dist.ncols();
    let nrow = dist.nrows();
    Mat::from_fn(nrow, ncol, |i, j| {
        let x = dist.get(i, j);
        f64::exp(-((x * epsilon).powi(2)))
    })
}

/// Bump Radial Basis function
/// Will set dist >= 1 / epsilon to 0, i.e., is a sparse RBF
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
pub fn rbf_bump_mat(dist: MatRef<'_, f64>, epsilon: &f64) -> Mat<f64> {
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
pub fn rbf_inverse_quadratic(dist: &[f64], epsilon: &f64) -> Vec<f64> {
    dist.par_iter()
        .map(|x| 1.0 / (1.0 + (*epsilon * x).powi(2)))
        .collect()
}

/// Inverse quadratic RBF for matrices
pub fn rbf_inverse_quadratic_mat(dist: MatRef<'_, f64>, epsilon: &f64) -> Mat<f64> {
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

/// Calculates the topological overlap measure for a given affinity matrix
/// Assumes a symmetric affinity matrix. Has the option to calculate the
/// signed and unsigned version. Supports both Version1 and Version2 algorithms.
pub fn calc_tom(affinity_mat: MatRef<'_, f64>, signed: bool, tom_type: TomType) -> Mat<f64> {
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

/// Function that assesses the cluster stability. Assumes that the rows
/// are the features and the columns the different resamples/bootstraps.
/// Returns a vector of tuples with the first being the average Jaccard
/// and the second the standard deviation.
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
