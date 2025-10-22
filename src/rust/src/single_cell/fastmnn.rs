#![allow(dead_code)]

use faer::{Mat, MatRef};
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};

use crate::single_cell::sc_knn_snn::*;

////////////
// Params //
////////////

/// Parameters for fastMNN batch correction
///
/// ### Fields
///
/// * `k` - Number of mutual nearest neighbours to identify
/// * `sigma` - Bandwidth of the Gaussian smoothing kernel (as proportion of
///   space radius after optional cosine normalisation)
/// * `knn_method` - Approximate nearest neighbour search method. One of
///   `"annoy"` or `"hnsw"`.
/// * `dist_metric` - Distance metric to use. One of `"cosine"` or
///   `"euclidean"`.
/// * `annoy_n_trees` - Number of trees for Annoy index (only used if
///   knn_method="annoy")
/// * `annoy_search_budget` - Search budget per tree for Annoy (only used if
///   knn_method="annoy")
/// * `cos_norm` - Apply cosine normalisation before computing distances
/// * `var_adj` - Apply variance adjustment to avoid kissing effects
/// * `svd_dim` - Number of dimensions for biological subspace. If 0,
///   biological variance correction is disabled.
#[derive(Clone, Debug)]
pub struct FastMnnParams {
    pub k: usize,
    pub sigma: f32,
    pub knn_method: String,
    pub dist_metric: String,
    pub annoy_n_trees: usize,
    pub annoy_search_budget: usize,
    pub cos_norm: bool,
    pub var_adj: bool,
    pub svd_dim: usize,
}

/////////////
// Helpers //
/////////////

/// Find mutual nearest neighbours from two KNN graphs
///
/// ### Params
///
/// * `left_knn` - KNN indices for left batch (each row is cell's neighbours)
/// * `right_knn` - KNN indices for right batch
///
/// ### Returns
///
/// (left_indices, right_indices) of MNN pairs
pub fn find_mutual_nns(
    left_knn: &[Vec<usize>],
    right_knn: &[Vec<usize>],
) -> (Vec<usize>, Vec<usize>) {
    let right_sets: Vec<FxHashSet<usize>> = right_knn
        .iter()
        .map(|neighbours| neighbours.iter().copied().collect())
        .collect();

    left_knn
        .par_iter()
        .enumerate()
        .fold(
            || (Vec::new(), Vec::new()),
            |(mut left_mnn, mut right_mnn), (left_idx, left_neighbours)| {
                for &right_idx in left_neighbours {
                    if right_sets[right_idx].contains(&left_idx) {
                        left_mnn.push(left_idx);
                        right_mnn.push(right_idx);
                    }
                }
                (left_mnn, right_mnn)
            },
        )
        .reduce(
            || (Vec::new(), Vec::new()),
            |(mut l1, mut r1), (l2, r2)| {
                l1.extend(l2);
                r1.extend(r2);
                (l1, r1)
            },
        )
}

/// Compute raw correction vectors from MNN pairs
///
/// ### Params
///
/// * `data_1` - Left batch data (cells x genes)
/// * `data_2` - Right batch data (cells x genes)
/// * `mnn_1` - MNN indices in left batch
/// * `mnn_2` - MNN indices in right batch
///
/// ### Returns
///
/// Matrix of correction vectors averaged per unique cell in right batch
/// (cells x genes)
pub fn compute_correction_vecs(
    data_1: &MatRef<f32>,
    data_2: &MatRef<f32>,
    mnn_1: &[usize],
    mnn_2: &[usize],
) -> Mat<f32> {
    let n_features = data_1.ncols();
    let ncells_2 = data_2.nrows();

    let mut accum: FxHashMap<usize, (Vec<f32>, usize)> = FxHashMap::default();

    for (&idx1, &idx2) in mnn_1.iter().zip(mnn_2.iter()) {
        let (sums, count) = accum
            .entry(idx2)
            .or_insert_with(|| (vec![0_f32; n_features], 0));
        for g in 0..n_features {
            sums[g] += data_1.get(idx1, g) - data_2.get(idx2, g);
        }
        *count += 1;
    }

    let mut averaged = Mat::zeros(ncells_2, n_features);
    for (cell_idx, (sums, count)) in accum.iter() {
        let n = *count as f32;
        for g in 0..n_features {
            averaged[(*cell_idx, g)] = sums[g] / n;
        }
    }

    averaged
}

/// Logspace addition to avoid underflow
#[inline]
fn logspace_add(log_a: f32, log_b: f32) -> f32 {
    if log_a.is_infinite() && log_a.is_sign_negative() {
        return log_b;
    }
    if log_b.is_infinite() && log_b.is_sign_negative() {
        return log_a;
    }

    let max = log_a.max(log_b);
    max + ((log_a - max).exp() + (log_b - max).exp()).ln()
}

/// Smooth correction vectors using Gaussian kernel weighted by MNN density
///
/// ### Params
///
/// * `averaged` - Averaged correction vectors (cells_with_mnn x features)
/// * `mnn_indices` - Indices of cells that have MNN pairs
/// * `data_2` - Data of second batch. (Cells x features.)
/// * `sigma_2` - Bandwidth squared
///
/// ### Returns
///
/// Smoothed correction vectors for all cells (cells x features)
pub fn smooth_gaussian_kernel_mnn(
    averaged: &MatRef<f32>,
    mnn_indices: &[usize],
    data_2: &MatRef<f32>,
    sigma_square: f32,
) -> Mat<f32> {
    let n_cells = data_2.nrows();
    let n_features_dist = data_2.ncols();
    let n_features = averaged.ncols();
    let inv_sigma_square = 1.0 / sigma_square;

    let (output, log_total_prob) = (0..mnn_indices.len())
        .into_par_iter()
        .fold(
            || {
                (
                    vec![vec![0_f32; n_features]; n_cells],
                    vec![f32::NEG_INFINITY; n_cells],
                )
            },
            |(mut output, mut log_total_prob), mnn_i| {
                let mnn_cell_idx = mnn_indices[mnn_i];
                let mut log_probs = vec![0_f32; n_cells];

                // Compute log weights
                for other in 0..n_cells {
                    let mut dist_2 = 0_f32;
                    for g in 0..n_features_dist {
                        let diff = data_2.get(mnn_cell_idx, g) - data_2.get(other, g);
                        dist_2 += diff * diff;
                    }
                    log_probs[other] = -dist_2 * inv_sigma_square;
                }

                // Compute density
                let density = mnn_indices
                    .iter()
                    .map(|&idx| log_probs[idx])
                    .fold(f32::NEG_INFINITY, logspace_add);

                // Update output
                for other in 0..n_cells {
                    let log_weight = log_probs[other] - density;
                    let weight = log_weight.exp();

                    log_total_prob[other] = logspace_add(log_total_prob[other], log_weight);

                    for g in 0..n_features {
                        output[other][g] += averaged.get(mnn_i, g) * weight;
                    }
                }
                (output, log_total_prob)
            },
        )
        .reduce(
            || {
                (
                    vec![vec![0_f32; n_features]; n_cells],
                    vec![f32::NEG_INFINITY; n_cells],
                )
            },
            |(mut o1, mut p1), (o2, p2)| {
                for i in 0..n_cells {
                    for g in 0..n_features {
                        o1[i][g] += o2[i][g];
                    }
                    p1[i] = logspace_add(p1[i], p2[i]);
                }
                (o1, p1)
            },
        );

    // Normalise
    let mut result = Mat::zeros(n_cells, n_features);
    for i in 0..n_cells {
        let norm = log_total_prob[i].exp();
        if norm > 0_f32 {
            for g in 0..n_features {
                result[(i, g)] = output[i][g] / norm;
            }
        }
    }
    result
}

/// Merge two batches
///
/// ### Params
///
/// * `data_1` - First batch (cells x features)
/// * `data_2` - Second batch (cells x features)
/// * `params` - `FastMnnParams` params with all of the parameters for this
///   run
/// * `seed` - Random seed for reproducibility
/// * `verbose` - Controls verbosity of the function
///
/// ### Returns
///
/// Corrected data_2 (genes x features)
pub fn merge_two_batches(
    data_1: &MatRef<f32>,
    data_2: &MatRef<f32>,
    params: &FastMnnParams,
    seed: usize,
    verbose: bool,
) -> Mat<f32> {
    let sigma_square = params.sigma * params.sigma;

    let knn_method: KnnSearch = parse_knn_method(&params.knn_method).unwrap_or(KnnSearch::Annoy);

    let (knn_1_to_2, knn_2_to_1) = match knn_method {
        KnnSearch::Hnsw => {
            let index_2 = build_hnsw_index(*data_2, &params.dist_metric, seed);
            let (knn_1_to_2, _) = query_hnsw_index(
                *data_1,
                &index_2,
                &params.dist_metric,
                params.k,
                false,
                verbose,
            );

            let index_1 = build_hnsw_index(*data_1, &params.dist_metric, seed);
            let (knn_2_to_1, _) = query_hnsw_index(
                *data_2,
                &index_1,
                &params.dist_metric,
                params.k,
                false,
                verbose,
            );

            (knn_1_to_2, knn_2_to_1)
        }
        KnnSearch::Annoy => {
            let index_2 = build_annoy_index(*data_2, params.annoy_n_trees, seed);
            let (knn_1_to_2, _) = query_annoy_index(
                *data_1,
                &index_2,
                &params.dist_metric,
                params.k,
                params.annoy_search_budget,
                false,
                verbose,
            );

            let index_1 = build_annoy_index(*data_1, params.annoy_n_trees, seed);
            let (knn_2_to_1, _) = query_annoy_index(
                *data_2,
                &index_1,
                &params.dist_metric,
                params.k,
                params.annoy_search_budget,
                false,
                verbose,
            );

            (knn_1_to_2, knn_2_to_1)
        }
    };

    let (mnn_1, mnn_2) = find_mutual_nns(&knn_1_to_2, &knn_2_to_1);

    if mnn_1.is_empty() {
        if verbose {
            eprintln!("Warning: No MNN pairs found");
        }
        return data_2.to_owned();
    }

    if verbose {
        println!("Found {} MNN pairs", mnn_1.len());
    }

    let averaged = compute_correction_vecs(data_1, data_2, &mnn_1, &mnn_2);
    let corrections = smooth_gaussian_kernel_mnn(&averaged.as_ref(), &mnn_2, data_2, sigma_square);

    let mut corrected = data_2.to_owned();
    for cell in 0..data_2.nrows() {
        for gene in 0..data_2.ncols() {
            corrected[(cell, gene)] += corrections[(cell, gene)];
        }
    }

    corrected
}

/// Fast MNN with cell order tracking
///
/// ### Params
///
/// * `batches` - Vec of PCA matrices per batch (cells x n_pcs)
/// * `original_indices` - Vec of original cell indices per batch
/// * `params` - `FastMnnParams` params with all of the parameters for this
///   run
/// * `seed` - Random seed for reproducibility
/// * `verbose` - Controls verbosity of the function
///
/// ### Returns
///
/// (corrected_pca, output_to_original_mapping)
pub fn fast_mnn(
    batches: Vec<Mat<f32>>,
    original_indices: Vec<Vec<usize>>,
    params: &FastMnnParams,
    seed: usize,
    verbose: bool,
) -> (Mat<f32>, Vec<usize>) {
    assert_eq!(batches.len(), original_indices.len());

    let mut merged = batches[0].to_owned();
    let mut index_map = original_indices[0].clone();

    for (batch, batch_indices) in batches
        .into_iter()
        .zip(original_indices.into_iter())
        .skip(1)
    {
        merged = merge_two_batches(&merged.as_ref(), &batch.as_ref(), params, seed, verbose);
        index_map.extend(batch_indices);
    }

    (merged.to_owned(), index_map)
}

/// Reorder corrected PCA back to original cell order
///
/// ### Params
///
/// * `corrected_pca` - Output from fast_mnn (cells x n_pcs)
/// * `output_to_original` - Mapping from output row -> original index
///
/// ### Returns
///
/// Reordered matrix matching original cell order
pub fn reorder_to_original(corrected_pca: &Mat<f32>, output_to_original: &[usize]) -> Mat<f32> {
    let n_cells = corrected_pca.nrows();
    let n_pcs = corrected_pca.ncols();

    // Create inverse mapping: original_idx -> output_idx
    let mut original_to_output = vec![0; n_cells];
    for (output_idx, &original_idx) in output_to_original.iter().enumerate() {
        original_to_output[original_idx] = output_idx;
    }

    // Reorder rows
    Mat::from_fn(n_cells, n_pcs, |row, col| {
        *corrected_pca.get(original_to_output[row], col)
    })
}
