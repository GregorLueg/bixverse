use faer::{MatRef, RowRef};
use instant_distance::{Builder, HnswMap, Point as DistancePoint, Search};
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::time::Instant;
use thousands::Separable;

use crate::core::data::sparse_structures::*;
use crate::core::graph::annoy::AnnoyIndex;
use crate::core::graph::knn::{parse_ann_dist, AnnDist};

///////////
// Enums //
///////////

/// SNN similarity method
#[derive(Clone, Copy)]
pub enum SnnSimilarityMethod {
    /// This will calculate the Jaccard similarity as weight
    Intersection,
    /// This will calculate the Rank version as a weight
    Rank,
}

/// Enum for the different methods
pub enum KnnSearch {
    /// Annoy-based
    Annoy,
    /// Hierarchical Navigable Small World
    Hnsw,
}

////////////////
// Structures //
////////////////

/// Point structure for the HNSW index
#[derive(Clone, Debug)]
pub struct Point(Vec<f32>, AnnDist);

impl DistancePoint for Point {
    /// Distance function. This is Euclidean distance without squaring for
    /// speed gains. Does not change the rank order in KNN generation.
    fn distance(&self, other: &Self) -> f32 {
        match self.1 {
            // No & needed, Copy does the work
            AnnDist::Euclidean => {
                let mut sum = 0.0f32;
                for i in 0..self.0.len() {
                    let diff = self.0[i] - other.0[i];
                    sum += diff * diff;
                }
                sum
            }
            AnnDist::Cosine => {
                let mut dot = 0.0f32;
                let mut norm_a = 0.0f32;
                let mut norm_b = 0.0f32;

                for i in 0..self.0.len() {
                    dot += self.0[i] * other.0[i];
                    norm_a += self.0[i] * self.0[i];
                    norm_b += other.0[i] * other.0[i];
                }

                1.0 - (dot / (norm_a.sqrt() * norm_b.sqrt()))
            }
        }
    }
}

/////////////
// Helpers //
/////////////

/// Compute distance between two points
///
/// Helper function to quickly calculate the implemented distances additionally
///
/// ### Params
///
/// * `a` - RowRef to cell a.
/// * `b` - RowRef to cell b.
///
/// ### Returns
///
/// The distance between the two cells based on the embedding.
pub fn compute_distance_knn(a: RowRef<f32>, b: RowRef<f32>, metric: &AnnDist) -> f32 {
    match metric {
        AnnDist::Euclidean => {
            let mut sum = 0.0f32;
            for i in 0..a.ncols() {
                let diff = a[i] - b[i];
                sum += diff * diff;
            }
            sum.sqrt()
        }
        AnnDist::Cosine => {
            let mut dot = 0.0f32;
            let mut norm_a = 0.0f32;
            let mut norm_b = 0.0f32;

            for i in 0..a.ncols() {
                dot += a[i] * b[i];
                norm_a += a[i] * a[i];
                norm_b += b[i] * b[i];
            }

            1.0 - (dot / (norm_a.sqrt() * norm_b.sqrt()))
        }
    }
}

/// Helper function to get the KNN method
///
/// ### Params
///
/// * `s` - Type of KNN algorithm to use
///
/// ### Returns
///
/// Option of the HvgMethod (some not yet implemented)
pub fn get_knn_method(s: &str) -> Option<KnnSearch> {
    match s.to_lowercase().as_str() {
        "annoy" => Some(KnnSearch::Annoy),
        "hnsw" => Some(KnnSearch::Hnsw),
        _ => None,
    }
}

/// Helper function to get the type of sNN similarity
///
/// ### Params
///
/// * `s` - Type of SNN similarity to use
///
/// ### Returns
///
/// Option of the SnnSimilarityMethod
pub fn get_snn_similiarity_method(s: &str) -> Option<SnnSimilarityMethod> {
    match s.to_lowercase().as_str() {
        "jaccard" => Some(SnnSimilarityMethod::Intersection),
        "rank" => Some(SnnSimilarityMethod::Rank),
        _ => None,
    }
}

/// Helper function to create a kNN mat with self
///
/// ### Params
///
/// * knn_graph - The kNN graph structure in which rows represent samples and
///   the columns represent the neighbours
///
/// ### Results
///
/// Updated version with self added
pub fn build_nn_map(knn_graph: &[Vec<usize>]) -> Vec<Vec<usize>> {
    (0..knn_graph.len())
        .map(|i| {
            let mut neighbors = knn_graph[i].clone();
            neighbors.push(i);
            neighbors
        })
        .collect()
}

/// Helper to calculate smooth kNN distances
///
/// ### Params
///
/// * `dist` - Slice of distance vectors from the kNN
/// * `k` - Number of neighbours.
/// * `local_connectivity` - Determines the minimum distance threshold (rho) by
///   interpolating between the nearest non-zero neighbors. If it's 1.5, you'd
///   interpolate between the 1st and 2nd neighbor distances.
/// * `smook_k_tol` - Tolerance parameter that controls:
///   - Whether to apply interpolation when computing rho
///   - When to stop the binary search (when |psum - target| < smook_k_tol)
/// * `min_k_dist_scale` - Minimum scaling factor for sigma to prevent it from
///   becoming too small
///
/// ### Returns
///
/// Tuple of (rho, sigma)
pub fn smooth_knn_dist(
    dist: &[Vec<f32>],
    k: f32,
    local_connectivity: f32,
    smook_k_tol: f32,
    min_k_dist_scale: f32,
) -> (Vec<f32>, Vec<f32>) {
    let n = dist.len();
    let n_neighbours = dist[0].len();
    let target = k.log2();

    let mean_dist = dist.iter().flat_map(|d| d.iter()).sum::<f32>() / (n * n_neighbours) as f32;

    let res: Vec<(f32, f32)> = dist
        .par_iter()
        .map(|dist_i| {
            let mut rho = 0.0_f32;

            let non_zero_dist: Vec<f32> = dist_i.iter().filter(|&&d| d > 0.0).copied().collect();

            if non_zero_dist.len() >= local_connectivity as usize {
                let index = local_connectivity.floor() as usize;
                let interpolation = local_connectivity - local_connectivity.floor();

                if index > 0 {
                    rho = non_zero_dist[index - 1];
                    if interpolation > smook_k_tol {
                        rho += interpolation * (non_zero_dist[index] - non_zero_dist[index - 1]);
                    }
                } else {
                    rho = interpolation * non_zero_dist[0];
                }
            } else if !non_zero_dist.is_empty() {
                rho = *non_zero_dist
                    .iter()
                    .max_by(|a, b| a.partial_cmp(b).unwrap())
                    .unwrap();
            }

            // binary search for sigma
            let mut lo = 0.0_f32;
            let mut hi = f32::INFINITY;
            let mut mid = 1.0_f32;

            for _ in 0..64 {
                let mut psum = 0.0_f32;

                for j in 1..n_neighbours {
                    let d = dist_i[j] - rho;
                    if d > 0.0_f32 {
                        psum += (-d / mid).exp();
                    } else {
                        psum += 1.0;
                    }
                }

                if (psum - target).abs() < smook_k_tol {
                    break;
                }

                if psum > target {
                    hi = mid;
                    mid = (lo + hi) / 2.0;
                } else {
                    lo = mid;
                    if hi == f32::INFINITY {
                        mid *= 2.0;
                    } else {
                        mid = (lo + hi) / 2.0;
                    }
                }
            }

            let mut sigma = mid;

            if rho > 0.0 {
                let mean_i = dist_i.iter().sum::<f32>() / n_neighbours as f32;
                if sigma < min_k_dist_scale * mean_i {
                    sigma = min_k_dist_scale * mean_i;
                }
            } else if sigma < min_k_dist_scale * mean_dist {
                sigma = min_k_dist_scale * mean_dist
            }

            (sigma, rho)
        })
        .collect();

    let sigmas = res.iter().map(|r| r.0).collect();
    let rhos = res.iter().map(|r| r.1).collect();

    (sigmas, rhos)
}

/// Helper function to transform kNN data into CompressedSparseData
///
/// ### Params
///
/// * `knn_indices` - The indices of the k-nearest neighbours.
/// * `knn_dists` - The distances to the k-nearest neighbours.
/// * `n_obs` - Number of observations in the data.
///
/// ### Return
///
/// `CompressedSparseData` in CSR format with distances to the k-nearest
/// neighbours stored.
pub fn knn_to_sparse_dist(
    knn_indices: &[Vec<usize>],
    knn_dists: &[Vec<f32>],
    n_obs: usize,
) -> CompressedSparseData<f32> {
    let mut rows = Vec::new();
    let mut cols = Vec::new();
    let mut vals = Vec::new();

    for i in 0..knn_indices.len() {
        for j in 0..knn_indices[i].len() {
            let neighbor = knn_indices[i][j];
            let dist = if neighbor == i { 0.0 } else { knn_dists[i][j] };

            if dist != 0.0 {
                rows.push(i);
                cols.push(neighbor);
                vals.push(dist);
            }
        }
    }

    coo_to_csr(&rows, &cols, &vals, (n_obs, n_obs))
}

///////////////////
// KNN functions //
///////////////////

/// Build an Annoy index
///
/// ### Params
///
/// * `mat` - The data matrix. Rows represent the samples, columns represent
///   the embedding dimensions
/// * `n_trees` - Number of trees to use to build the index
/// * `usize` - Random seed for reproducibility
///
/// ### Return
///
/// The `AnnoyIndex`.
pub fn build_annoy_index(mat: MatRef<f32>, n_trees: usize, seed: usize) -> AnnoyIndex {
    AnnoyIndex::new(mat, n_trees, seed)
}

/// Helper function to query a given Annoy index
///
/// ### Params
///
/// * `query_mat` - The query matrix containing the samples x features
/// * `dist_metric` - The distance metric to use. One of `"euclidean"` or
///   `"cosine"`.
/// * `k` - Number of neighbours to return
/// * `search_budget` - Search budget per tree
/// * `return_dist` - Shall the distances between the different points be
///   returned
///
/// ### Returns
///
/// A tuple of `(knn_indices, optional distances)`
pub fn query_annoy_index(
    query_mat: MatRef<f32>,
    index: &AnnoyIndex,
    dist_metric: &str,
    k: usize,
    search_budget: usize,
    return_dist: bool,
    verbose: bool,
) -> (Vec<Vec<usize>>, Option<Vec<Vec<f32>>>) {
    let n_samples = query_mat.nrows();
    let ann_dist = parse_ann_dist(dist_metric).unwrap();
    let search_k = Some(k * search_budget);
    let counter = Arc::new(AtomicUsize::new(0));

    if return_dist {
        let results: Vec<(Vec<usize>, Vec<f32>)> = (0..n_samples)
            .into_par_iter()
            .map(|i| {
                let (neighbors, dists) = index.query_row(query_mat.row(i), &ann_dist, k, search_k);
                if verbose {
                    let count = counter.fetch_add(1, Ordering::Relaxed) + 1;
                    if count % 100_000 == 0 {
                        println!(" Processed {} / {} cells.", count, n_samples);
                    }
                }
                (neighbors, dists)
            })
            .collect();
        let (indices, distances) = results.into_iter().unzip();
        (indices, Some(distances))
    } else {
        let indices: Vec<Vec<usize>> = (0..n_samples)
            .into_par_iter()
            .map(|i| {
                let (neighbors, _) = index.query_row(query_mat.row(i), &ann_dist, k, search_k);
                if verbose {
                    let count = counter.fetch_add(1, Ordering::Relaxed) + 1;
                    if count % 100_000 == 0 {
                        println!(" Processed {} / {} cells.", count, n_samples);
                    }
                }
                neighbors
            })
            .collect();
        (indices, None)
    }
}

/// Build HNSW index
///
/// ### Params
///
/// * `mat` - The data matrix. Rows represent the samples, columns represent
///   the embedding dimensions
/// * `dist_metric` - Which distance metric to use. One of `"euclidean"` or
///   `"cosine"`
/// * `usize` - Random seed for reproducibility
///
/// ### Return
///
/// The
pub fn build_hnsw_index(mat: MatRef<f32>, dist_metric: &str, seed: usize) -> HnswMap<Point, usize> {
    let ann_dist = parse_ann_dist(dist_metric).unwrap();
    let n_samples = mat.nrows();
    let points: Vec<Point> = (0..n_samples)
        .into_par_iter()
        .map(|i| Point(mat.row(i).iter().cloned().collect(), ann_dist))
        .collect();

    Builder::default()
        .seed(seed as u64)
        .build(points, (0..n_samples).collect::<Vec<_>>())
}

/// Query HNSW index
///
/// ### Params
///
/// * `query_mat` - The query matrix containing the samples x features
/// * `dist_metric` - The distance metric to use. One of `"euclidean"` or
///   `"cosine"`.
/// * `k` - Number of neighbours to return
/// * `return_dist` - Shall the distances between the different points be
///   returned
///
/// ### Returns
///
/// A tuple of `(knn_indices, optional distances)`
pub fn query_hnsw_index(
    query_mat: MatRef<f32>,
    index: &HnswMap<Point, usize>,
    dist_metric: &str,
    k: usize,
    return_dist: bool,
    verbose: bool,
) -> (Vec<Vec<usize>>, Option<Vec<Vec<f32>>>) {
    let n_samples = query_mat.nrows();
    let ann_dist = parse_ann_dist(dist_metric).unwrap();
    let n_query = query_mat.nrows();
    let counter = Arc::new(AtomicUsize::new(0));

    let results: Vec<_> = (0..n_query)
        .into_par_iter()
        .map(|i| {
            let point = Point(query_mat.row(i).iter().cloned().collect(), ann_dist);
            let mut search = Search::default();
            let neighbors: Vec<usize> = index
                .search(&point, &mut search)
                .take(k)
                .map(|item| *item.value)
                .collect();

            let dists = if return_dist {
                let mut dists = Vec::with_capacity(k);
                for &neighbor_idx in &neighbors {
                    let dist = compute_distance_knn(
                        query_mat.row(i),
                        query_mat.row(neighbor_idx),
                        &ann_dist,
                    );
                    dists.push(dist);
                }
                Some(dists)
            } else {
                None
            };

            let count = counter.fetch_add(1, Ordering::Relaxed) + 1;
            if verbose && count % 100_000 == 0 {
                println!(
                    " Processed {} / {} cells.",
                    count.separate_with_underscores(),
                    n_samples.separate_with_underscores()
                );
            }

            (neighbors, dists)
        })
        .collect();

    let indices = results.iter().map(|(i, _)| i.clone()).collect();
    let distances = if return_dist {
        Some(results.iter().map(|(_, d)| d.clone().unwrap()).collect())
    } else {
        None
    };

    (indices, distances)
}

////////////////////
// Main functions //
////////////////////

/// Get the kNN graph based on HNSW
///
/// This function generates the kNN graph via an approximate nearest neighbour
/// search based on the HNSW algorithm (hierarchical navigable small world).
///
/// ### Params
///
/// * `mat` - Matrix in which rows represent the samples and columns the
///   respective embeddings for that sample
/// * `no_neighbours` - Number of neighbours for the KNN graph
/// * `seed` - Seed for the HNSW algorithm
///
/// ### Returns
///
/// The k-nearest neighbours based on the HNSW algorithm
pub fn generate_knn_hnsw(
    mat: MatRef<f32>,
    dist_metric: &str,
    no_neighbours: usize,
    seed: usize,
    verbose: bool,
) -> Vec<Vec<usize>> {
    let start_index = Instant::now();

    let index = build_hnsw_index(mat, dist_metric, seed);

    let end_index = start_index.elapsed();
    if verbose {
        println!("Generated HNSW index: {:.2?}", end_index);
    }

    let start_search = Instant::now();

    let (indices, _) =
        query_hnsw_index(mat, &index, dist_metric, no_neighbours + 1, false, verbose);

    let res: Vec<Vec<usize>> = indices
        .into_iter()
        .enumerate()
        .map(|(i, mut neighbors)| {
            neighbors.retain(|&x| x != i);
            neighbors.truncate(no_neighbours);
            neighbors
        })
        .collect();

    let end_search = start_search.elapsed();
    if verbose {
        println!(
            "Identified approximate nearest neighbours via HNSW: {:.2?}.",
            end_search
        );
    }

    res
}

/// Get the kNN graph based on Annoy
///
/// This function generates the kNN graph based via an approximate nearest
/// neighbour search based on the Annoy algorithm (or a version thereof).
///
/// ### Params
///
/// * `mat` - Matrix in which rows represent the samples and columns the
///   respective embeddings for that sample
/// * `no_neighbours` - Number of neighbours for the KNN graph.
/// * `n_trees` - Number of trees to use for the search.
/// * `search_budget` - Search budget per given query.
/// * `seed` - Seed for the HNSW algorithm
///
/// ### Returns
///
/// The k-nearest neighbours based on the HNSW algorithm
pub fn generate_knn_annoy(
    mat: MatRef<f32>,
    dist_metric: &str,
    no_neighbours: usize,
    n_trees: usize,
    search_budget: usize,
    seed: usize,
    verbose: bool,
) -> Vec<Vec<usize>> {
    let start_index = Instant::now();

    let index = build_annoy_index(mat, n_trees, seed);

    let end_index = start_index.elapsed();
    if verbose {
        println!("Generated Annoy index: {:.2?}", end_index);
    }

    let start_search = Instant::now();

    let (indices, _) = query_annoy_index(
        mat,
        &index,
        dist_metric,
        no_neighbours + 1,
        search_budget,
        false,
        verbose,
    );

    let res: Vec<Vec<usize>> = indices
        .into_iter()
        .enumerate()
        .map(|(i, mut neighbors)| {
            neighbors.retain(|&x| x != i);
            neighbors.truncate(no_neighbours);
            neighbors
        })
        .collect();

    let end_search = start_search.elapsed();
    if verbose {
        println!(
            "Identified approximate nearest neighbours via Annoy: {:.2?}.",
            end_search
        );
    }

    res
}

/// Generate an sNN graph based on the kNN graph (full)
///
/// This version will compare all cells against all cells and generate an edge
/// if any neighbours are shared. This yields way denser graphs and is the
/// approach taken in the `bluster` R package to generate the sNN.
///
/// ### Params
///
/// * `knn_graph` - K-nearest neighbours data as a flat vector in column-major.
/// * `no_neighbours` - Number of neighbours in the kNN graph
/// * `pruning` - Below which Jaccard similarity to prune the edge. In this case
///   the weight is set to `0`.
/// * `method` - Which similarity method to use
/// * `verbose` - Controls verbosity of the function.
///
/// ### Returns
///
/// A tuple with `(<edges>, <weights>)`. The edges are stored in a way that the
/// the first edge points goes from the first element to the second, the second
/// edge from the third to the fourth, etc.
pub fn generate_snn_full(
    flat_knn: &[usize],
    k: usize,
    n_samples: usize,
    pruning: f32,
    method: SnnSimilarityMethod,
    verbose: bool,
) -> (Vec<usize>, Vec<f32>) {
    let mut reverse_mappings: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n_samples];

    let start_time = Instant::now();

    for i in 0..n_samples {
        reverse_mappings[i].push((i, 0));

        for neighbor_idx in 0..k {
            let neighbor = flat_knn[neighbor_idx * n_samples + i];
            reverse_mappings[neighbor].push((i, neighbor_idx + 1));
        }
    }

    let results: Vec<(usize, usize, f32)> = (0..n_samples)
        .into_par_iter()
        .flat_map(|j| {
            let mut scores = vec![0.0f32; n_samples];
            let mut added = Vec::new();

            for i in 0..=k {
                let cur_neighbor = if i == 0 {
                    j
                } else {
                    flat_knn[(i - 1) * n_samples + j]
                };

                for &(othernode, other_rank) in &reverse_mappings[cur_neighbor] {
                    if othernode < j {
                        match method {
                            SnnSimilarityMethod::Rank => {
                                let combined_rank = (i + other_rank) as f32;
                                if scores[othernode] == 0.0 {
                                    scores[othernode] = combined_rank;
                                    added.push(othernode);
                                } else if combined_rank < scores[othernode] {
                                    scores[othernode] = combined_rank;
                                }
                            }
                            SnnSimilarityMethod::Intersection => {
                                if scores[othernode] == 0.0 {
                                    added.push(othernode);
                                }
                                scores[othernode] += 1.0;
                            }
                        }
                    }
                }
            }

            added
                .into_iter()
                .filter_map(|othernode| {
                    let weight = match method {
                        SnnSimilarityMethod::Rank => {
                            let preliminary = k as f32 - scores[othernode] / 2.0;
                            let raw_weight = preliminary.max(1e-6);
                            raw_weight / k as f32
                        }
                        SnnSimilarityMethod::Intersection => {
                            scores[othernode] / (2.0 * (k as f32 + 1.0) - scores[othernode])
                        }
                    };

                    if weight >= pruning {
                        Some((j, othernode, weight))
                    } else {
                        None
                    }
                })
                .collect::<Vec<_>>()
        })
        .collect();

    let mut edges = Vec::with_capacity(results.len() * 2);
    let mut weights = Vec::with_capacity(results.len());

    for (i, j, weight) in results {
        edges.push(i);
        edges.push(j);
        weights.push(weight);
    }

    let end_snn = start_time.elapsed();

    if verbose {
        println!("Transformed kNN into a full sNN graph: {:.2?}", end_snn);
    }

    (edges, weights)
}

/// Generate an sNN graph based on the kNN graph (limited)
///
/// This version will only compare cells to the neighbouring cells and
/// deduplicate edges in taking the maximum weight between two given cells.
///
/// ### Params
///
/// * `knn_graph` - K-nearest neighbours data as a flat vector in column-major.
/// * `no_neighbours` - Number of neighbours in the kNN graph
/// * `pruning` - Below which Jaccard similarity to prune the edge. In this case
///   the weight is set to `0`.
/// * `method` - Which similarity method to use.
/// * `verbose` - Controls verbosity of the function.
///
/// ### Returns
///
/// A tuple with `(<edges>, <weights>)`. The edges are stored in a way that the
/// the first edge points goes from the first element to the second, the second
/// edge from the third to the fourth, etc.
pub fn generate_snn_limited(
    flat_knn: &[usize],
    k: usize,
    n_samples: usize,
    pruning: f32,
    method: SnnSimilarityMethod,
    verbose: bool,
) -> (Vec<usize>, Vec<f32>) {
    // We need to use a hashmap to store unique edges (smaller_idx, larger_idx)
    // -> weight
    let start_time = Instant::now();

    let edge_map: FxHashMap<(usize, usize), f32> = (0..n_samples)
        .into_par_iter()
        .flat_map(|i| {
            let mut edges = Vec::new();

            // only consider edges to this cell's k nearest neighbors
            for neighbor_idx in 0..k {
                let j = flat_knn[neighbor_idx * n_samples + i];

                // Calculate sNN similarity between cell i and its neighbor j
                let weight = match method {
                    SnnSimilarityMethod::Intersection => {
                        // Get neighbors of both cells
                        let neighbors_i: FxHashSet<usize> = (0..k)
                            .map(|idx| flat_knn[idx * n_samples + i])
                            .chain(std::iter::once(i)) // include self
                            .collect();

                        let neighbors_j: FxHashSet<usize> = (0..k)
                            .map(|idx| flat_knn[idx * n_samples + j])
                            .chain(std::iter::once(j)) // include self
                            .collect();

                        let intersection_count =
                            neighbors_i.intersection(&neighbors_j).count() as f32;
                        intersection_count / (2.0 * (k as f32 + 1.0) - intersection_count)
                        // Jaccard
                    }
                    SnnSimilarityMethod::Rank => {
                        // build ranks i
                        let mut ranks_i = FxHashMap::default();
                        ranks_i.insert(i, 0); // self at rank 0
                        for (rank, neighbor) in
                            (0..k).map(|idx| flat_knn[idx * n_samples + i]).enumerate()
                        {
                            ranks_i.insert(neighbor, rank + 1);
                        }

                        // build ranks j
                        let mut ranks_j = FxHashMap::default();
                        ranks_j.insert(j, 0); // self at rank 0
                        for (rank, neighbor) in
                            (0..k).map(|idx| flat_knn[idx * n_samples + j]).enumerate()
                        {
                            ranks_j.insert(neighbor, rank + 1);
                        }

                        // find minimum combined rank of shared neighbors
                        let min_combined_rank = ranks_i
                            .keys()
                            .filter(|&neighbor| ranks_j.contains_key(neighbor))
                            .map(|neighbor| ranks_i[neighbor] + ranks_j[neighbor])
                            .min()
                            .unwrap_or(2 * k)
                            as f32;

                        let preliminary = k as f32 - min_combined_rank / 2.0;
                        let raw_weight = preliminary.max(1e-6);

                        raw_weight / k as f32
                    }
                };

                if weight >= pruning {
                    // Store edge with smaller index first to ensure uniqueness
                    let edge_key = if i < j { (i, j) } else { (j, i) };
                    edges.push((edge_key, weight));
                }
            }

            edges
        })
        .collect::<Vec<_>>()
        .into_iter()
        .fold(FxHashMap::default(), |mut acc, (edge_key, weight)| {
            // Keep the maximum weight if we see the same edge multiple times
            acc.entry(edge_key)
                .and_modify(|existing_weight| {
                    if weight > *existing_weight {
                        *existing_weight = weight;
                    }
                })
                .or_insert(weight);
            acc
        });

    let mut edges = Vec::with_capacity(edge_map.len() * 2);
    let mut weights = Vec::with_capacity(edge_map.len());

    for ((i, j), weight) in edge_map {
        edges.push(i);
        edges.push(j);
        weights.push(weight);
    }

    let end_snn = start_time.elapsed();

    if verbose {
        println!("Transformed kNN into an sNN graph: {:.2?}", end_snn);
    }

    (edges, weights)
}
