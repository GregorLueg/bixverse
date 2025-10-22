#![allow(dead_code)]

use extendr_api::List;
use faer::{Mat, MatRef};
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};

use crate::single_cell::processing::*;
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
/// * `no_pcs` - Number of PCs to use for the MNN calculations
/// * `random_svd` - Boolean. Shall randomised SVD be used.
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
    pub no_pcs: usize,
    pub random_svd: bool,
}

impl FastMnnParams {
    /// Generate the FastMnnParams from an R list
    ///
    /// Should values not be found within the List, the parameters will default
    /// to sensible defaults.
    ///
    /// ### Params
    ///
    /// * `r_list` - The list with the fastMNN parameters.
    ///
    /// ### Return
    ///
    /// The `FastMnnParams` with all of the parameters.
    pub fn from_r_list(r_list: List) -> Self {
        let fastmnn_list = r_list.into_hashmap();

        let k = fastmnn_list
            .get("k")
            .and_then(|v| v.as_integer())
            .unwrap_or(20) as usize;

        let sigma = fastmnn_list
            .get("sigma")
            .and_then(|v| v.as_real())
            .unwrap_or(0.1) as f32;

        let knn_method = std::string::String::from(
            fastmnn_list
                .get("knn_method")
                .and_then(|v| v.as_str())
                .unwrap_or("annoy"),
        );

        let dist_metric = std::string::String::from(
            fastmnn_list
                .get("dist_metric")
                .and_then(|v| v.as_str())
                .unwrap_or("cosine"),
        );

        let annoy_n_trees = fastmnn_list
            .get("annoy_n_trees")
            .and_then(|v| v.as_integer())
            .unwrap_or(100) as usize;

        let annoy_search_budget = fastmnn_list
            .get("annoy_search_budget")
            .and_then(|v| v.as_integer())
            .unwrap_or(100) as usize;

        let cos_norm = fastmnn_list
            .get("cos_norm")
            .and_then(|v| v.as_logical())
            .map(|rb| rb.is_true())
            .unwrap_or(true);

        let var_adj = fastmnn_list
            .get("var_adj")
            .and_then(|v| v.as_logical())
            .map(|rb| rb.is_true())
            .unwrap_or(true);

        let no_pcs = fastmnn_list
            .get("no_pcs")
            .and_then(|v| v.as_integer())
            .unwrap_or(30) as usize;

        let random_svd = fastmnn_list
            .get("random_svd")
            .and_then(|v| v.as_logical())
            .map(|rb| rb.is_true())
            .unwrap_or(true);

        Self {
            k,
            sigma,
            knn_method,
            dist_metric,
            annoy_n_trees,
            annoy_search_budget,
            cos_norm,
            var_adj,
            no_pcs,
            random_svd,
        }
    }
}

/////////////
// Helpers //
/////////////

/// Apply cosine normalisation (L2 normalisation) to each row
///
/// ### Params
///
/// * `mat` - The matrix on which to apply the Cosine normalisation per row
///
/// ### Returns
///
/// Per row L2-normalised data
pub fn cosine_normalise(mat: &Mat<f32>) -> Mat<f32> {
    Mat::from_fn(mat.nrows(), mat.ncols(), |row, col| {
        let norm = mat.get(row, ..).norm_l2();
        if norm > 1e-8 {
            mat[(row, col)] / norm
        } else {
            0.0
        }
    })
}

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

    // Validate indices
    for &idx in mnn_indices {
        assert!(
            idx < n_cells,
            "mnn_indices contains {} but data_2 only has {} rows",
            idx,
            n_cells
        );
    }

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
                        output[other][g] += averaged.get(mnn_cell_idx, g) * weight;
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
    let mut corrections =
        smooth_gaussian_kernel_mnn(&averaged.as_ref(), &mnn_2, data_2, sigma_square);

    // In merge_two_batches, after computing corrections:
    if params.var_adj {
        let scaling = adjust_shift_variance(
            data_1,
            data_2,
            &corrections.as_ref(),
            sigma_square,
            None, // or pass restricted indices
            None,
        );

        for cell in 0..corrections.nrows() {
            for gene in 0..corrections.ncols() {
                corrections[(cell, gene)] *= scaling[cell].max(1.0);
            }
        }
    }

    let mut corrected = data_2.to_owned();
    for cell in 0..data_2.nrows() {
        for gene in 0..data_2.ncols() {
            corrected[(cell, gene)] += corrections[(cell, gene)];
        }
    }

    // CHANGE: Concatenate data_1 and corrected data_2
    let n_cells_total = data_1.nrows() + corrected.nrows();
    let n_features = data_1.ncols();

    Mat::from_fn(n_cells_total, n_features, |row, col| {
        if row < data_1.nrows() {
            data_1[(row, col)]
        } else {
            corrected[(row - data_1.nrows(), col)]
        }
    })
}

/// Compute squared distance from point to line defined by ref + t*grad
fn sq_distance_to_line_buf(ref_point: &[f32], grad: &[f32], point: &[f32]) -> f32 {
    let scale: f32 = ref_point
        .iter()
        .zip(grad.iter())
        .zip(point.iter())
        .map(|((r, g), p)| (r - p) * g)
        .sum();

    ref_point
        .iter()
        .zip(grad.iter())
        .zip(point.iter())
        .map(|((r, g), p)| {
            let diff = r - p;
            let perp = diff - scale * g;
            perp * perp
        })
        .sum()
}

/// Adjust shift variance to avoid kissing effects
///
/// ### Params
///
/// * `data1` - Reference batch (cells x genes)
/// * `data2` - Target batch (cells x genes)
/// * `corrections` - Correction vectors (cells x genes)
/// * `sigma2` - Bandwidth squared
/// * `restrict1` - Optional subset of cells from batch1 (None = all)
/// * `restrict2` - Optional subset of cells from batch2 (None = all)
///
/// ### Returns
///
/// Scaling factors for each cell in batch2
pub fn adjust_shift_variance(
    data1: &MatRef<f32>,
    data2: &MatRef<f32>,
    corrections: &MatRef<f32>,
    sigma2: f32,
    restrict1: Option<&[usize]>,
    restrict2: Option<&[usize]>,
) -> Vec<f32> {
    let ngenes = data1.ncols();
    let ncells2 = data2.nrows();

    assert_eq!(data1.ncols(), ngenes);
    assert_eq!(data2.ncols(), ngenes);
    assert_eq!(corrections.nrows(), ncells2);
    assert_eq!(corrections.ncols(), ngenes);

    let restrict1_indices: Vec<usize> = restrict1
        .map(|r| r.to_vec())
        .unwrap_or_else(|| (0..data1.nrows()).collect());
    let restrict2_indices: Vec<usize> = restrict2
        .map(|r| r.to_vec())
        .unwrap_or_else(|| (0..ncells2).collect());

    let inv_sigma2 = 1.0 / sigma2;
    let mut output = vec![1.0f32; ncells2];

    // Pre-allocate working buffers - reused across iterations
    let mut grad = vec![0.0f32; ngenes];
    let mut curcell = vec![0.0f32; ngenes];
    let mut othercell = vec![0.0f32; ngenes];
    let mut distance1: Vec<(f32, f32)> = vec![(0.0, 0.0); restrict1_indices.len()];

    for cell in 0..ncells2 {
        // Copy correction vector
        for g in 0..ngenes {
            grad[g] = corrections[(cell, g)];
        }

        let l2norm: f32 = grad.iter().map(|&x| x * x).sum::<f32>().sqrt();

        if l2norm < 1e-8 {
            continue;
        }

        let inv_l2norm = 1.0 / l2norm;
        for g in grad.iter_mut() {
            *g *= inv_l2norm;
        }

        // Copy current cell
        for g in 0..ngenes {
            curcell[g] = data2[(cell, g)];
        }

        let curproj: f32 = grad.iter().zip(curcell.iter()).map(|(g, c)| g * c).sum();

        // Batch2 cumulative probability
        let mut prob2 = f32::NEG_INFINITY;
        let mut totalprob2 = f32::NEG_INFINITY;

        for &same_idx in &restrict2_indices {
            let (log_prob, should_add) = if same_idx == cell {
                (0.0, true)
            } else {
                for g in 0..ngenes {
                    othercell[g] = data2[(same_idx, g)];
                }
                let samedist = sq_distance_to_line_buf(&curcell, &grad, &othercell);
                let sameproj: f32 = grad.iter().zip(othercell.iter()).map(|(g, c)| g * c).sum();
                (-samedist * inv_sigma2, sameproj <= curproj)
            };

            totalprob2 = logspace_add(totalprob2, log_prob);
            if should_add {
                prob2 = logspace_add(prob2, log_prob);
            }
        }

        prob2 -= totalprob2;

        // Batch1 distances
        for (idx, &other_idx) in restrict1_indices.iter().enumerate() {
            for g in 0..ngenes {
                othercell[g] = data1[(other_idx, g)];
            }

            let proj: f32 = grad.iter().zip(othercell.iter()).map(|(g, c)| g * c).sum();
            let dist = sq_distance_to_line_buf(&curcell, &grad, &othercell);
            distance1[idx] = (proj, -dist * inv_sigma2);
        }

        if distance1.is_empty() {
            continue;
        }

        distance1.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        let totalprob1 = distance1
            .iter()
            .map(|(_, lw)| *lw)
            .fold(f32::NEG_INFINITY, logspace_add);

        let target = prob2 + totalprob1;
        let mut cumulative = f32::NEG_INFINITY;
        let mut ref_quan = distance1.last().unwrap().0;

        for &(proj, log_weight) in &distance1 {
            cumulative = logspace_add(cumulative, log_weight);
            if cumulative >= target {
                ref_quan = proj;
                break;
            }
        }

        output[cell] = (ref_quan - curproj) * inv_l2norm;
    }

    output
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
    let n_pcs = corrected_pca.ncols();
    let n_output_cells = corrected_pca.nrows();

    let n_original_cells = output_to_original.iter().copied().max().unwrap_or(0) + 1;

    assert_eq!(
        n_output_cells,
        output_to_original.len(),
        "Mismatch: corrected_pca has {} rows but output_to_original has {} elements",
        n_output_cells,
        output_to_original.len()
    );

    let mut original_to_output = vec![0; n_original_cells];
    for (output_idx, &original_idx) in output_to_original.iter().enumerate() {
        original_to_output[original_idx] = output_idx;
    }

    Mat::from_fn(n_original_cells, n_pcs, |row, col| {
        *corrected_pca.get(original_to_output[row], col)
    })
}

/// Split PCA matrix by batch indices
///
/// ### Params
///
/// * `pca_all` - Full PCA matrix (all cells x n_pcs)
/// * `batch_indices` - Batch assignment for each cell
///
/// ### Returns
///
/// (matrices per batch, original cell indices per batch)
pub fn split_pca_by_batch(
    pca_all: &Mat<f32>,
    batch_indices: &[usize],
) -> (Vec<Mat<f32>>, Vec<Vec<usize>>) {
    let n_features = pca_all.ncols();
    let n_batches = batch_indices.iter().copied().max().unwrap_or(0) + 1;

    let mut batch_cells: Vec<Vec<usize>> = vec![Vec::new(); n_batches];
    for (cell_idx, &batch) in batch_indices.iter().enumerate() {
        batch_cells[batch].push(cell_idx);
    }

    let batches: Vec<Mat<f32>> = batch_cells
        .iter()
        .map(|cells| {
            Mat::from_fn(cells.len(), n_features, |row, col| {
                pca_all[(cells[row], col)]
            })
        })
        .collect();

    (batches, batch_cells)
}

///////////////////
// Main function //
///////////////////

/// Perform fastMNN batch correction on single cell data
///
/// Computes PCA on all cells, splits by batch, applies fastMNN correction,
/// and returns the corrected PCA matrix in the original cell order.
///
/// ### Params
///
/// * `f_path` - Path to single cell data file
/// * `cell_indices` - Indices of cells to include
/// * `gene_indices` - Indices of genes to include
/// * `batch_indices` - Batch assignment for each cell
/// * `params` - FastMNN parameters
/// * `verbose` - Controls verbosity of the function
/// * `seed` - Random seed for reproducibility
///
/// ### Returns
///
/// Batch-corrected PCA matrix (cells x n_pcs) in original cell order
pub fn fast_mnn_main(
    f_path: &str,
    cell_indices: &[usize],
    gene_indices: &[usize],
    batch_indices: &[usize],
    params: &FastMnnParams,
    verbose: bool,
    seed: usize,
) -> Mat<f32> {
    // calculate the PCA across everything
    let (pca_all, _, _, _) = pca_on_sc(
        f_path,
        cell_indices,
        gene_indices,
        params.no_pcs,
        params.random_svd,
        seed,
        false,
        verbose,
    );

    let (mut pca_batches, original_indices) = split_pca_by_batch(&pca_all, batch_indices);

    if params.cos_norm {
        for batch in pca_batches.iter_mut() {
            *batch = cosine_normalise(batch);
        }
    }

    let (corrected, index_map) = fast_mnn(pca_batches, original_indices, params, seed, verbose);

    reorder_to_original(&corrected, &index_map)
}
