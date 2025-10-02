use faer::{Mat, MatRef};
use rayon::prelude::*;
use std::collections::BinaryHeap;

use crate::assert_symmetric_mat;
use crate::core::base::cors_similarity::*;
use crate::core::base::utils::*;
use crate::core::graph::knn::SimilarityItem;

/////////////
// Helpers //
/////////////

/// Create an affinity matrix from a distance matrix
///
/// Applies a scaled exponential similarity kernel based on K-nearest neighbors.
/// The kernel uses an adaptive sigma parameter calculated from the average
/// distance to K nearest neighbors. Sounds fancy, doesn't it... ?
///
/// ### Params
///
/// * `dist` - The distance matrix.
/// * `k` - Number of neighbours to consider
/// * `mu` - Controls the Gaussian kernel strength, i.e., sigma parameter
///
/// ### Returns
///
/// The resulting affinity matrix.
fn affinity_from_distance(dist: &MatRef<f64>, k: usize, mu: f64) -> Result<Mat<f64>, String> {
    let n = dist.nrows();

    // compute average distance to K nearest neighbors for each sample (parallelized)
    let knn_avg_dist: Vec<f64> = (0..n)
        .into_par_iter()
        .map(|i| {
            let mut distances: Vec<f64> = (0..n)
                .filter(|&j| i != j)
                .map(|j| *dist.get(i, j))
                .collect();

            distances.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());

            let k_actual = k.min(distances.len());
            distances[..k_actual].iter().sum::<f64>() / k_actual as f64
        })
        .collect();

    // apply gaussian kernel with scaled sigma
    let mut affinity = Mat::zeros(n, n);

    for i in 0..n {
        affinity[(i, i)] = 1.0;
        for j in (i + 1)..n {
            let dij = *dist.get(i, j);
            let sigma = mu * (knn_avg_dist[i] + knn_avg_dist[j] + dij) / 3.0;

            if sigma < f64::EPSILON {
                affinity[(i, j)] = 0.0;
                affinity[(j, i)] = 0.0;
            } else {
                let weight = (-dij.powi(2) / (2.0 * sigma.powi(2))).exp();
                affinity[(i, j)] = weight;
                affinity[(j, i)] = weight;
            }
        }
    }

    Ok(affinity)
}

/// Normalise the data across rows
///
/// ### Params
///
/// * `mat` - The matrix to normalise.
///
/// ### Returns
///
/// The (row-)normalised matrix.
fn normalise_rows(mat: &MatRef<f64>) -> Mat<f64> {
    let (nrows, ncols) = mat.shape();
    Mat::from_fn(nrows, ncols, |i, j| {
        let row_sum: f64 = (0..ncols).map(|k| mat.get(i, k)).sum();
        if row_sum > 0.0 {
            mat.get(i, j) / row_sum
        } else {
            0.0
        }
    })
}

/// KNN thresholding for the matrix
///
/// For each sample, keeps only the K most similar neighbors and normalises
/// their weights to sum to 1. Used to emphasise local structure.
///
/// ### Params
///
/// * `mat` - The similarity matrix to threshold
/// * `k` - Number of nearest neighbors to retain per sample
///
/// ### Returns
///
/// A sparse matrix where each row contains at most K non-zero entries.
fn knn_threshold(mat: &MatRef<f64>, k: usize) -> Mat<f64> {
    let n = mat.nrows();
    let mut result = Mat::zeros(n, n);

    let rows: Vec<Vec<(usize, f64)>> = (0..n)
        .into_par_iter()
        .map(|i| {
            let mut heap = BinaryHeap::with_capacity(k + 1);

            for j in 0..n {
                if i != j {
                    let sim = *mat.get(i, j);
                    heap.push(SimilarityItem {
                        index: j,
                        similarity: sim,
                    });
                    if heap.len() > k {
                        heap.pop();
                    }
                }
            }

            let neighbors: Vec<_> = heap
                .into_iter()
                .map(|item| (item.index, item.similarity))
                .collect();

            let knn_sum: f64 = neighbors.iter().map(|(_, v)| v).sum();

            neighbors
                .into_iter()
                .map(|(idx, val)| (idx, val / knn_sum))
                .collect()
        })
        .collect();

    for (i, neighbors) in rows.iter().enumerate() {
        for &(j, normalized_val) in neighbors {
            result[(i, j)] = normalized_val;
        }
    }

    result
}

/// B0 normalisation for SNF update step
///
/// Applies modified normalisation where diagonal elements are set to 0.5 and
/// off-diagonal elements are scaled by alpha. Increases stability during the
/// fusion process
///
/// ### Params
///
/// * `mat` - The matrix to normalise.
/// * `alpha` - Normalisation factor for off-diagonal elements (typically 1.0)
///
/// ### Returns
///
/// The normalised matrix with diagonal set to 0.5
fn b0_normalise(mat: &MatRef<f64>, alpha: f64) -> Mat<f64> {
    let n = mat.nrows();
    let normalised = normalise_rows(mat);

    Mat::from_fn(n, n, |i, j| {
        if i == j {
            0.5
        } else {
            normalised.get(i, j) / (2.0 * alpha)
        }
    })
}

/// Calculate Gower distance between columns
///
/// I need this one, as I cannot just transpose the other one.
///
/// ### Params
///
/// * `mat` - The data matrix (features × samples)
/// * `is_cat` - Boolean vector indicating which rows (features) are categorical
/// * `ranges` - Optional pre-computed ranges for continuous variables
///
/// ### Returns
///
/// The Gower distance matrix with values in [0, 1]
pub fn snf_gower_dist(
    mat: &MatRef<f64>,
    is_cat: &[bool],
    ranges: Option<&[f64]>,
) -> Result<Mat<f64>, String> {
    let (nrows, ncols) = mat.shape();

    if is_cat.len() != nrows {
        return Err(format!(
            "is_cat length {} doesn't match features {}",
            is_cat.len(),
            nrows
        ));
    }

    // compute ranges in parallel
    let computed_ranges: Vec<f64> = if let Some(r) = ranges {
        if r.len() != nrows {
            return Err(format!(
                "ranges length {} doesn't match features {}",
                r.len(),
                nrows
            ));
        }
        r.to_vec()
    } else {
        (0..nrows)
            .into_par_iter()
            .map(|i| {
                if is_cat[i] {
                    1.0
                } else {
                    let mut min_val = f64::INFINITY;
                    let mut max_val = f64::NEG_INFINITY;
                    for j in 0..ncols {
                        let val = *mat.get(i, j);
                        min_val = min_val.min(val);
                        max_val = max_val.max(val);
                    }
                    let range = max_val - min_val;
                    if range < f64::EPSILON {
                        1.0
                    } else {
                        range
                    }
                }
            })
            .collect()
    };

    let mut res = Mat::zeros(ncols, ncols);

    let pairs: Vec<(usize, usize)> = (0..ncols)
        .flat_map(|i| ((i + 1)..ncols).map(move |j| (i, j)))
        .collect();

    let results: Vec<(usize, usize, f64)> = pairs
        .par_iter()
        .map(|&(i, j)| {
            let mut total_dist = 0.0;

            for k in 0..nrows {
                let val_i = *mat.get(k, i);
                let val_j = *mat.get(k, j);

                let dist = if is_cat[k] {
                    if (val_i - val_j).abs() < f64::EPSILON {
                        0.0
                    } else {
                        1.0
                    }
                } else {
                    (val_i - val_j).abs() / computed_ranges[k]
                };

                total_dist += dist;
            }

            (i, j, total_dist / nrows as f64)
        })
        .collect();

    for (i, j, dist) in results {
        res[(i, j)] = dist;
        res[(j, i)] = dist;
    }

    Ok(res)
}

////////////////////
// Main functions //
////////////////////

/// Generates an affinity matrix for SNF on continuous values
///
/// ### Params
///
/// * `data` - The underlying data. Assumes the orientation features x samples!
/// * `distance_type` - One of the implemented distances.
/// * `k` - Number of neighbours to consider
/// * `mu` - Controls the Gaussian kernel strength
/// * `normalise` - Shall the data be normalised prior to distance calculation
///
/// ### Returns
///
/// The affinity matrix based on continuous values
pub fn make_affinity_continuous(
    data: &MatRef<f64>,
    distance_type: &str,
    k: usize,
    mu: f64,
    normalise: bool,
) -> Result<Mat<f64>, String> {
    let dist_type = parse_distance_type(distance_type)
        .ok_or_else(|| format!("Invalid distance type: {}", distance_type))?;

    let normalised_data = if normalise {
        scale_matrix_col(data, true)
    } else {
        data.to_owned()
    };

    let dist_mat = match dist_type {
        DistanceType::L1Norm => column_pairwise_l1_norm(&normalised_data.as_ref()),
        DistanceType::L2Norm => column_pairwise_l2_norm(&normalised_data.as_ref()),
        DistanceType::Cosine => column_pairwise_cosine_dist(&normalised_data.as_ref()),
        DistanceType::Canberra => column_pairwise_canberra_dist(&normalised_data.as_ref()),
    };

    affinity_from_distance(&dist_mat.as_ref(), k, mu)
}

/// Generates an affinity matrix for SNF on mixed feature types
///
/// ### Params
///
/// * `data` - The underlying data. Assumes the orientation features x samples!
/// * `is_cat` - Which of the features are categorical.
/// * `k` - Number of neighbours to consider
/// * `mu` - Controls the Gaussian kernel strength
///
/// ### Returns
///
/// The affinity matrix based on Gower distance
pub fn make_affinity_mixed(
    data: &MatRef<f64>,
    is_cat: &[bool],
    k: usize,
    mu: f64,
) -> Result<Mat<f64>, String> {
    let dist_mat = snf_gower_dist(data, is_cat, None)?;
    affinity_from_distance(&dist_mat.as_ref(), k, mu)
}

/// Generates an affinity matrix for SNF on categorical features
///
/// ### Params
///
/// * `data` - The underlying data. Assumes the orientation features x samples!
/// * `k` - Number of neighbours to consider
/// * `mu` - Controls the Gaussian kernel strength
///
/// ### Returns
///
/// The affinity matrix based on Gower distance
pub fn make_affinity_categorical(
    data: &MatRef<i32>,
    k: usize,
    mu: f64,
) -> Result<Mat<f64>, String> {
    let dist_mat = column_pairwise_hamming_cat(data);
    affinity_from_distance(&dist_mat.as_ref(), k, mu)
}

/// Run similarity network fusion
///
/// Add more details @Claude
///
/// ### Params
///
/// * `aff_mats` - Slice of matrix references representing the individual
///   affinity matrices. The dimensions need to be the same and they need to
///   be symmetric.
/// * `k` - Number of neighbours to consider.
/// * `t` - Number of iterations to run the algorithm for.
/// * `alpha` - Normalisation hyperparameter controlling fusion strength
///   (typically 1.0)
///
/// ### Returns
///
/// Final results
pub fn snf(aff_mats: &[MatRef<f64>], k: usize, t: usize, alpha: f64) -> Result<Mat<f64>, String> {
    if aff_mats.is_empty() {
        return Err("At least one affinity matrix required".to_string());
    }

    let n = aff_mats[0].nrows();

    for (i, mat) in aff_mats.iter().enumerate() {
        assert_symmetric_mat!(mat);
        if mat.nrows() != n {
            return Err(format!(
                "Matrix {} has different size: {} x {} (expected {} x {})",
                i,
                mat.nrows(),
                mat.ncols(),
                n,
                n
            ));
        }
    }

    let m = aff_mats.len();

    let mut aff = aff_mats
        .par_iter()
        .map(normalise_rows)
        .collect::<Vec<Mat<f64>>>();

    let wk = aff
        .par_iter()
        .map(|mat_i| knn_threshold(&mat_i.as_ref(), k))
        .collect::<Vec<Mat<f64>>>();

    for _ in 0..t {
        for v in 0..m {
            // compute sum of all P matrices except the current one
            let mut p_sum: Mat<f64> = Mat::zeros(n, n);
            for (idx, mat) in aff.iter().enumerate() {
                if idx != v {
                    p_sum += mat;
                }
            }
            p_sum = &p_sum / (m - 1) as f64;

            let fused = &wk[v] * &p_sum * wk[v].transpose();
            aff[v] = b0_normalise(&fused.as_ref(), alpha)
        }
    }

    // final fusion - mean across all matrices
    let mut w: Mat<f64> = aff
        .par_iter()
        .cloned()
        .reduce(|| Mat::zeros(n, n), |acc, mat| acc + mat);

    w = &w / m as f64;

    w = normalise_rows(&w.as_ref());

    w = (&w + w.transpose()) / 2.0;

    for i in 0..n {
        w[(i, i)] = 1.0;
    }

    Ok(w)
}
