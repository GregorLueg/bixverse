use extendr_api::List;
use faer::MatRef;
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use statrs::distribution::ChiSquared;
use statrs::distribution::ContinuousCDF;

use crate::core::data::sparse_structures::*;
use crate::core::graph::knn::{parse_ann_dist, AnnDist};
use crate::single_cell::sc_knn_snn::*;
use crate::utils::general::MatSliceView;

// use crate::core::base::stats::hypergeom_pval;

//////////
// kBET //
//////////

/// Calculate kBET-based mixing scores on kNN data
///
/// ### Params
///
/// * `knn_data` - KNN data. Outer vector represents the cells, while the inner
///   vector represents
/// * `batches` - Vector indicating the batches.
///
/// ### Return
///
/// Numerical vector indicating with the p-values from the ChiSquare test
pub fn kbet(knn_data: &Vec<Vec<usize>>, batches: &Vec<usize>) -> Vec<f64> {
    let mut batch_counts = FxHashMap::default();
    for &batch in batches {
        *batch_counts.entry(batch).or_insert(0) += 1;
    }
    let total = batches.len() as f64;
    let batch_ids: Vec<usize> = batch_counts.keys().copied().collect();
    let dof = (batch_ids.len() - 1) as f64;

    knn_data
        .par_iter()
        .map(|neighbours| {
            let k = neighbours.len() as f64;
            let mut neighbours_count = FxHashMap::default();
            for &neighbour_idx in neighbours {
                *neighbours_count.entry(batches[neighbour_idx]).or_insert(0) += 1;
            }

            // Chi-square test: Σ (observed - expected)² / expected
            let mut chi_square = 0.0;
            for &batch_id in &batch_ids {
                let expected = k * (batch_counts[&batch_id] as f64 / total);
                let observed = *neighbours_count.get(&batch_id).unwrap_or(&0) as f64;
                chi_square += (observed - expected).powi(2) / expected;
            }

            // Compute p-value from chi-square distribution
            1.0 - ChiSquared::new(dof).unwrap().cdf(chi_square)
        })
        .collect()
}

///////////
// BBKNN //
///////////

// Constants needed for the UMAP stuff
const SMOOTH_K_TOLERANCE: f32 = 1e-5;
const MIN_K_DIST_SCALE: f32 = 1e-3;

#[derive(Clone, Debug)]
pub struct BbknnParams {
    pub neighbours_within_batch: usize,
    pub knn_method: String,
    pub dist_metric: String,
    pub set_op_mix_ratio: f32,
    pub local_connectivity: f32,
    pub annoy_n_trees: usize,
    pub annoy_search_budget: usize,
    pub trim: Option<usize>,
}

impl BbknnParams {
    pub fn from_r_list(r_list: List) -> Self {
        let bbknn_list = r_list.into_hashmap();

        let neighbours_within_batch = bbknn_list
            .get("neighbours_within_batch")
            .and_then(|v| v.as_integer())
            .unwrap_or(3) as usize;

        let knn_method = std::string::String::from(
            bbknn_list
                .get("knn_method")
                .and_then(|v| v.as_str())
                .unwrap_or("annoy"),
        );

        let dist_metric = std::string::String::from(
            bbknn_list
                .get("dist_metric")
                .and_then(|v| v.as_str())
                .unwrap_or("cosine"),
        );

        let set_op_mix_ratio = bbknn_list
            .get("set_op_mix_ratio")
            .and_then(|v| v.as_real())
            .unwrap_or(1.0) as f32;

        let local_connectivity = bbknn_list
            .get("local_connectivity")
            .and_then(|v| v.as_real())
            .unwrap_or(1.0) as f32;

        let annoy_n_trees = bbknn_list
            .get("annoy_n_trees")
            .and_then(|v| v.as_integer())
            .unwrap_or(100) as usize;

        let annoy_search_budget = bbknn_list
            .get("annoy_search_budget")
            .and_then(|v| v.as_integer())
            .unwrap_or(100) as usize;

        let trim = bbknn_list
            .get("trim")
            .and_then(|v| v.as_integer())
            .unwrap_or(10 * neighbours_within_batch as i32) as usize;

        Self {
            neighbours_within_batch,
            knn_method,
            dist_metric,
            set_op_mix_ratio,
            local_connectivity,
            annoy_n_trees,
            annoy_search_budget,
            trim: Some(trim),
        }
    }
}

///////////////////
// Helpers BBKNN //
///////////////////

/// The function pulls out the batch specific indices of the cells
///
/// ### Params
///
/// ### Returns
fn get_batch_balanced_knn(
    mat: MatRef<f32>,
    batch_labels: &[usize], // Changed to usize
    knn_method: &KnnSearch,
    bbknn_params: &BbknnParams,
    seed: usize,
    verbose: bool,
) -> (Vec<Vec<usize>>, Vec<Vec<f32>>) {
    let n_cells = mat.nrows();

    // get unique batches
    let unique_batches: Vec<usize> = {
        let mut batches: Vec<_> = batch_labels.to_vec();
        batches.sort_unstable();
        batches.dedup();
        batches
    };

    let dist_metric: AnnDist = parse_ann_dist(&bbknn_params.dist_metric).unwrap_or(AnnDist::Cosine);

    let n_batches = unique_batches.len();

    let mut all_indices = vec![vec![0; bbknn_params.neighbours_within_batch * n_batches]; n_cells];
    let mut all_distances =
        vec![vec![0.0; bbknn_params.neighbours_within_batch * n_batches]; n_cells];
    let col_indices: Vec<usize> = (0..mat.ncols()).collect();

    for (batch_idx, &batch) in unique_batches.iter().enumerate() {
        if verbose {
            println!(
                "Processing batch {} / {}: {}",
                batch_idx + 1,
                n_batches,
                batch
            );
        }

        let batch_cell_indices: Vec<usize> = batch_labels
            .iter()
            .enumerate()
            .filter(|(_, b)| **b == batch)
            .map(|(i, _)| i)
            .collect();

        let sub_matrix = MatSliceView::new(mat, &batch_cell_indices, &col_indices);
        let sub_matrix = sub_matrix.to_owned();

        let (neighbor_indices, neighbour_dists) = match knn_method {
            KnnSearch::Annoy => {
                // annoy path with updated functions
                let index =
                    build_annoy_index(sub_matrix.as_ref(), bbknn_params.annoy_n_trees, seed);
                query_annoy_index(
                    mat,
                    &index,
                    &bbknn_params.dist_metric,
                    bbknn_params.neighbours_within_batch + 1,
                    bbknn_params.annoy_search_budget,
                    false,
                    verbose,
                )
            }
            KnnSearch::Hnsw => {
                // hnsw path with updated functions
                let index = build_hnsw_index(sub_matrix.as_ref(), &bbknn_params.dist_metric, seed);
                query_hnsw_index(
                    mat,
                    &index,
                    &bbknn_params.dist_metric,
                    bbknn_params.neighbours_within_batch + 1,
                    false,
                    verbose,
                )
            }
        };

        // PRINT HERE - after querying, before the loop
        println!(
            "Batch {}: neighbor_indices[0] = {:?}",
            batch_idx, &neighbor_indices[0]
        );
        println!(
            "Batch {}: batch_cell_indices = {:?}",
            batch_idx,
            &batch_cell_indices[..10.min(batch_cell_indices.len())]
        );

        let col_start = batch_idx * bbknn_params.neighbours_within_batch;

        for cell_idx in 0..n_cells {
            let mut added = 0;
            let mut k_idx = 0;

            while added < bbknn_params.neighbours_within_batch {
                let local_idx = neighbor_indices[cell_idx][k_idx];
                let global_idx = batch_cell_indices[local_idx];

                // Skip self-loops
                if global_idx != cell_idx {
                    let dist = match dist_metric {
                        AnnDist::Euclidean => mat
                            .row(cell_idx)
                            .iter()
                            .zip(mat.row(global_idx).iter())
                            .map(|(a, b)| (a - b).powi(2))
                            .sum::<f32>()
                            .sqrt(),
                        AnnDist::Cosine => {
                            1.0 - mat
                                .row(cell_idx)
                                .iter()
                                .zip(mat.row(global_idx).iter())
                                .map(|(a, b)| a * b)
                                .sum::<f32>()
                                / (mat.row(cell_idx).iter().map(|x| x * x).sum::<f32>().sqrt()
                                    * mat
                                        .row(global_idx)
                                        .iter()
                                        .map(|x| x * x)
                                        .sum::<f32>()
                                        .sqrt())
                        }
                    };

                    all_indices[cell_idx][col_start + added] = global_idx;
                    all_distances[cell_idx][col_start + added] = dist;
                    added += 1;
                }

                k_idx += 1;
            }
        }

        // PRINT HERE - after the population loop
        println!(
            "After batch {}, all_indices[0] = {:?}",
            batch_idx, &all_indices[0]
        );
        println!(
            "After batch {}, all_distances[0] = {:?}",
            batch_idx, &all_distances[0]
        );
    }

    (all_indices, all_distances)
}

/// Sorts the kNN by distance
fn sort_knn_by_distance(knn_indices: &mut [Vec<usize>], knn_dists: &mut [Vec<f32>]) {
    for i in 0..knn_indices.len() {
        let mut pairs: Vec<_> = knn_dists[i]
            .iter()
            .zip(knn_indices[i].iter())
            .map(|(d, idx)| (*d, *idx))
            .collect();

        pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

        for (j, (dist, idx)) in pairs.into_iter().enumerate() {
            knn_dists[i][j] = dist;
            knn_indices[i][j] = idx;
        }
    }
}

/// Compute membership strengths
fn compute_membership_strengths(
    knn_indices: &[Vec<usize>],
    knn_dists: &[Vec<f32>],
    sigmas: &[f32],
    rhos: &[f32],
) -> (Vec<usize>, Vec<usize>, Vec<f32>) {
    let n_samples = knn_indices.len();
    let n_neighbours = knn_indices[0].len();

    let mut rows = Vec::with_capacity(n_samples * n_neighbours);
    let mut cols = Vec::with_capacity(n_samples * n_neighbours);
    let mut vals = Vec::with_capacity(n_samples * n_neighbours);

    for i in 0..n_samples {
        for j in 0..n_neighbours {
            let neighbor = knn_indices[i][j];

            let val = if neighbor == i {
                0.0
            } else if knn_dists[i][j] - rhos[i] <= 0.0 || sigmas[i] == 0.0 {
                1.0
            } else {
                (-(knn_dists[i][j] - rhos[i]) / sigmas[i]).exp()
            };

            rows.push(i);
            cols.push(neighbor);
            vals.push(val);
        }
    }

    (rows, cols, vals)
}

fn apply_set_operations(
    connectivities: CompressedSparseData<f32>,
    set_op_mix_ratio: f32,
) -> CompressedSparseData<f32> {
    let (nrow, ncol) = connectivities.shape;
    let transpose = {
        let csc_transpose = connectivities.transform(); // A^T in CSC
                                                        // Manually reinterpret CSC of A^T as CSR
        let mut coo_rows = Vec::new();
        let mut coo_cols = Vec::new();
        let mut coo_vals = Vec::new();

        for col in 0..csc_transpose.indptr.len() - 1 {
            for idx in csc_transpose.indptr[col]..csc_transpose.indptr[col + 1] {
                let row = csc_transpose.indices[idx];
                // For A^T: swap row/col to get proper CSR representation
                coo_rows.push(col);
                coo_cols.push(row);
                coo_vals.push(csc_transpose.data[idx]);
            }
        }
        coo_to_csr(&coo_rows, &coo_cols, &coo_vals, (ncol, nrow))
    };

    // Element-wise multiply: A .* A^T
    let prod = sparse_multiply_elementwise(&connectivities, &transpose);

    // set_op_mix_ratio * (A + A^T - A.*A^T) + (1 - set_op_mix_ratio) * (A.*A^T)
    let union_part = sparse_add_csr(&connectivities, &transpose);
    let union_part = sparse_subtract_csr(&union_part, &prod);
    let union_part = sparse_scalar_multiply_csr(&union_part, set_op_mix_ratio);

    let intersect_part = sparse_scalar_multiply_csr(&prod, 1.0 - set_op_mix_ratio);

    let res = sparse_add_csr(&union_part, &intersect_part);

    eliminate_zeros(res)
}

fn trim_graph(
    mut connectivities: CompressedSparseData<f32>,
    trim: usize,
) -> CompressedSparseData<f32> {
    let n = connectivities.shape.0;
    let mut thresholds = vec![0.0f32; n];

    // Compute thresholds
    for i in 0..n {
        let row_start = connectivities.indptr[i];
        let row_end = connectivities.indptr[i + 1];
        let row_data = &connectivities.data[row_start..row_end];

        if row_data.len() <= trim {
            continue;
        }

        let mut sorted = row_data.to_vec();
        sorted.sort_by(|a, b| b.partial_cmp(a).unwrap());
        thresholds[i] = sorted[trim - 1];
    }

    // Apply trimming twice (row then column)
    for _ in 0..2 {
        for i in 0..n {
            let row_start = connectivities.indptr[i];
            let row_end = connectivities.indptr[i + 1];

            for j in row_start..row_end {
                if connectivities.data[j] < thresholds[i] {
                    connectivities.data[j] = 0.0;
                }
            }
        }

        connectivities = eliminate_zeros(connectivities);
        connectivities = connectivities.transform();
    }

    connectivities
}

////////////////
// Main BBKNN //
////////////////

pub fn bbknn(
    mat: MatRef<f32>,
    batch_labels: &[usize],
    bbknn_params: &BbknnParams,
    seed: usize,
    verbose: bool,
) -> (CompressedSparseData<f32>, CompressedSparseData<f32>) {
    // parse it and worst case, I default to Annoy
    let knn_method = get_knn_method(&bbknn_params.knn_method).unwrap_or(KnnSearch::Annoy);

    // 1. Get batch-balanced k-NN
    let (mut knn_indices, mut knn_dists) =
        get_batch_balanced_knn(mat, batch_labels, &knn_method, bbknn_params, seed, verbose);

    // 2. Sort the distance by KNN <- Revisit if needed...
    sort_knn_by_distance(&mut knn_indices, &mut knn_dists);

    // After sort_knn_by_distance
    println!("After sorting, knn_indices[0] = {:?}", &knn_indices[0]);
    println!("After sorting, knn_dists[0] = {:?}", &knn_dists[0]);

    // 3. Compute UMAP connectivities
    let n_neighbours = knn_indices[0].len();
    let (sigmas, rhos) = smooth_knn_dist(
        &knn_dists,
        n_neighbours as f32,
        bbknn_params.local_connectivity,
        SMOOTH_K_TOLERANCE,
        MIN_K_DIST_SCALE,
    );

    let (rows, cols, vals) = compute_membership_strengths(&knn_indices, &knn_dists, &sigmas, &rhos);
    let n_obs = mat.nrows();

    let mut connectivities = coo_to_csr(&rows, &cols, &vals, (n_obs, n_obs));

    // 4. Apply set operations
    connectivities = apply_set_operations(connectivities, bbknn_params.set_op_mix_ratio);

    // 5. Create the distance matrix
    let dist = knn_to_sparse_dist(&knn_indices, &knn_dists, n_obs);

    // 6. Trimming
    if let Some(trim_val) = bbknn_params.trim {
        if trim_val > 0 {
            connectivities = trim_graph(connectivities, trim_val);
        }
    }

    (dist, connectivities)
}
