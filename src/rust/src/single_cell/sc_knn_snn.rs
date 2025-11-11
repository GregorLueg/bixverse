use extendr_api::List;
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
use crate::core::graph::nn_descent::NNDescent;
use crate::single_cell::cell_aggregations::{MetaCellParams, SuperCellParams};
use crate::single_cell::methods::seacells::SEACellsParams;
use crate::single_cell::methods::vision_hotspot::HotSpotParams;

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
    /// NNDescent
    NNDescent,
}

////////////
// Params //
////////////

/// KnnParams
///
/// ### Fields
///
/// **General**
///
/// * `knn_method` - Which of the kNN methods to use. One of `"annoy"`, `"hnsw"`
///   or `"nndescent"`.
/// * `ann_dist` - Approximate nearest neighbour distance measure. One of
///   `"euclidean"` or `"cosine"`.
/// * `k` - Number of neighbours to search
///
/// **Annoy**
///
/// * `n_tree` - Number of trees for the generation of the index
/// * `search_budget` - Search budget during querying
///
/// **NN Descent**
///
/// * `max_iter` - Maximum iterations for the algorithm
/// * `rho` - Sampling rate for the algorithm
/// * `delta` - Early termination criterium
pub struct KnnParams {
    // general params
    pub knn_method: String,
    pub ann_dist: String,
    pub k: usize,
    // annoy params
    pub n_tree: usize,
    pub search_budget: usize,
    // nn descent params
    pub max_iter: usize,
    pub rho: f32,
    pub delta: f32,
}

impl KnnParams {
    /// Generate KnnParams from an R list
    ///
    /// Should values not be found within the List, the parameters will default
    /// to sensible defaults based on heuristics.
    ///
    /// ### Params
    ///
    /// * `r_list` - The list with the Boost parameters.
    ///
    /// ### Returns
    ///
    /// The `KnnParams` with all parameters set.
    pub fn from_r_list(r_list: List) -> Self {
        let params_list = r_list.into_hashmap();

        // general
        let knn_method = std::string::String::from(
            params_list
                .get("knn_method")
                .and_then(|v| v.as_str())
                .unwrap_or("annoy"),
        );

        let ann_dist = std::string::String::from(
            params_list
                .get("ann_dist")
                .and_then(|v| v.as_str())
                .unwrap_or("cosine"),
        );

        let k = params_list
            .get("k")
            .and_then(|v| v.as_integer())
            .unwrap_or(15) as usize;

        // annoy
        let n_tree = params_list
            .get("n_tree")
            .and_then(|v| v.as_integer())
            .unwrap_or(100) as usize;

        let search_budget = params_list
            .get("search_budget")
            .and_then(|v| v.as_integer())
            .unwrap_or(100) as usize;

        // nn descent
        let max_iter = params_list
            .get("nn_max_iter")
            .and_then(|v| v.as_integer())
            .unwrap_or(25) as usize;

        let rho = params_list
            .get("rho")
            .and_then(|v| v.as_real())
            .unwrap_or(1.0) as f32;

        let delta = params_list
            .get("delta")
            .and_then(|v| v.as_real())
            .unwrap_or(0.001) as f32;

        Self {
            knn_method,
            ann_dist,
            k,
            n_tree,
            search_budget,
            max_iter,
            rho,
            delta,
        }
    }

    /// Extract the kNN params from the SEACells params
    ///
    /// ### Params
    ///
    /// * `seacell_params` - Structure holding the `SEACellsParams`.
    ///
    /// ### Returns
    ///
    /// The `KnnParams` with all parameters set.
    pub fn from_sea_cells_params(seacell_params: &SEACellsParams) -> Self {
        Self {
            knn_method: seacell_params.knn_method.clone(),
            ann_dist: seacell_params.ann_dist.clone(),
            k: seacell_params.k,
            n_tree: seacell_params.n_trees,
            search_budget: seacell_params.search_budget,
            max_iter: seacell_params.nn_max_iter,
            rho: seacell_params.rho,
            delta: seacell_params.delta,
        }
    }

    /// Extract the kNN params from the MetaCell parameters.
    ///
    /// ### Params
    ///
    /// * `metacell_params` - Structure holding the `MetaCellParams`.
    ///
    /// ### Returns
    ///
    /// The `KnnParams` with all parameters set.
    pub fn from_metacell_params(metacell_params: &MetaCellParams) -> Self {
        Self {
            knn_method: metacell_params.knn_method.clone(),
            ann_dist: metacell_params.ann_dist.clone(),
            k: metacell_params.k,
            n_tree: metacell_params.n_trees,
            search_budget: metacell_params.search_budget,
            max_iter: metacell_params.nn_max_iter,
            rho: metacell_params.rho,
            delta: metacell_params.delta,
        }
    }

    /// Extract the kNN params from the SuperCells parameters.
    ///
    /// ### Params
    ///
    /// * `supercell_params` - Structure holding the `SuperCellParams`.
    ///
    /// ### Returns
    ///
    /// The `KnnParams` with all parameters set.
    pub fn from_supercell_params(supercell_params: &SuperCellParams) -> Self {
        Self {
            knn_method: supercell_params.knn_method.clone(),
            ann_dist: supercell_params.ann_dist.clone(),
            k: supercell_params.k,
            n_tree: supercell_params.n_trees,
            search_budget: supercell_params.search_budget,
            max_iter: supercell_params.nn_max_iter,
            rho: supercell_params.rho,
            delta: supercell_params.delta,
        }
    }

    /// Extract the kNN params from the HotSpot parameters.
    ///
    /// ### Params
    ///
    /// * `hotspot_params` - Structure holding the `HotSpotParams`.
    ///
    /// ### Returns
    ///
    /// The `KnnParams` with all parameters set.
    pub fn from_hotspot_params(hotspot_params: &HotSpotParams) -> Self {
        Self {
            knn_method: hotspot_params.knn_method.clone(),
            ann_dist: hotspot_params.ann_dist.clone(),
            k: hotspot_params.k,
            n_tree: hotspot_params.n_tree,
            search_budget: hotspot_params.search_budget,
            max_iter: hotspot_params.max_iter,
            rho: hotspot_params.rho,
            delta: hotspot_params.delta,
        }
    }
}

////////////////
// Structures //
////////////////

/// Point structure for the HNSW index
#[derive(Clone, Debug)]
pub struct Point(Vec<f32>, AnnDist);

impl DistancePoint for Point {
    /// Distance function.
    ///
    /// ### Params
    ///
    /// * `other` - The other point to compare to
    ///
    /// ### Returns
    ///
    /// The distance between self and the other point.
    #[inline(always)]
    fn distance(&self, other: &Self) -> f32 {
        debug_assert_eq!(self.0.len(), other.0.len());

        // moaaaar unsafe and raw pointers...
        unsafe {
            let len = self.0.len();
            let ptr_a = self.0.as_ptr();
            let ptr_b = other.0.as_ptr();

            match self.1 {
                AnnDist::Euclidean => {
                    let mut sum = 0.0f32;
                    for i in 0..len {
                        let diff = *ptr_a.add(i) - *ptr_b.add(i);
                        sum += diff * diff;
                    }
                    sum
                }
                AnnDist::Cosine => {
                    let mut dot = 0.0f32;
                    let mut norm_a = 0.0f32;
                    let mut norm_b = 0.0f32;

                    for i in 0..len {
                        let a = *ptr_a.add(i);
                        let b = *ptr_b.add(i);
                        dot += a * b;
                        norm_a += a * a;
                        norm_b += b * b;
                    }

                    1.0 - (dot / (norm_a.sqrt() * norm_b.sqrt()))
                }
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
#[inline(always)]
pub fn compute_distance_knn(a: RowRef<f32>, b: RowRef<f32>, metric: &AnnDist) -> f32 {
    let ncols = a.ncols();

    // fast, unsafe path for contiguous memory
    if a.col_stride() == 1 && b.col_stride() == 1 {
        unsafe {
            let a_ptr = a.as_ptr();
            let b_ptr = b.as_ptr();

            match metric {
                AnnDist::Euclidean => {
                    let mut sum = 0.0f32;
                    for i in 0..ncols {
                        let diff = *a_ptr.add(i) - *b_ptr.add(i);
                        sum += diff * diff;
                    }
                    sum.sqrt()
                }
                AnnDist::Cosine => {
                    let mut dot = 0.0f32;
                    let mut norm_a = 0.0f32;
                    let mut norm_b = 0.0f32;

                    for i in 0..ncols {
                        let av = *a_ptr.add(i);
                        let bv = *b_ptr.add(i);
                        dot += av * bv;
                        norm_a += av * av;
                        norm_b += bv * bv;
                    }

                    1.0 - (dot / (norm_a.sqrt() * norm_b.sqrt()))
                }
            }
        }
    } else {
        // fallback
        match metric {
            AnnDist::Euclidean => {
                let mut sum = 0.0f32;
                for i in 0..ncols {
                    let diff = a[i] - b[i];
                    sum += diff * diff;
                }
                sum.sqrt()
            }
            AnnDist::Cosine => {
                let mut dot = 0.0f32;
                let mut norm_a = 0.0f32;
                let mut norm_b = 0.0f32;

                for i in 0..ncols {
                    dot += a[i] * b[i];
                    norm_a += a[i] * a[i];
                    norm_b += b[i] * b[i];
                }

                1.0 - (dot / (norm_a.sqrt() * norm_b.sqrt()))
            }
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
pub fn parse_knn_method(s: &str) -> Option<KnnSearch> {
    match s.to_lowercase().as_str() {
        "annoy" => Some(KnnSearch::Annoy),
        "hnsw" => Some(KnnSearch::Hnsw),
        "nndescent" => Some(KnnSearch::NNDescent),
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

///////////
// Annoy //
///////////

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
                    if count.is_multiple_of(100_000) {
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
                    if count.is_multiple_of(100_000) {
                        println!(" Processed {} / {} cells.", count, n_samples);
                    }
                }
                neighbors
            })
            .collect();
        (indices, None)
    }
}

//////////
// HNSW //
//////////

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
    let counter = Arc::new(AtomicUsize::new(0));

    let results: Vec<_> = (0..n_samples)
        .into_par_iter()
        .map(|i| {
            let point = Point(query_mat.row(i).iter().cloned().collect(), ann_dist);
            let mut search = Search::default();

            // capture both indices and distances in one pass
            let search_results: Vec<_> = index.search(&point, &mut search).take(k).collect();

            let neighbors: Vec<usize> = search_results.iter().map(|item| *item.value).collect();

            let dists = if return_dist {
                Some(search_results.iter().map(|item| item.distance).collect())
            } else {
                None
            };

            if verbose {
                let count = counter.fetch_add(1, Ordering::Relaxed) + 1;
                if count.is_multiple_of(100_000) {
                    println!(
                        " Processed {} / {} cells.",
                        count.separate_with_underscores(),
                        n_samples.separate_with_underscores()
                    );
                }
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

///////////////
// NNDescent //
///////////////

/// Get the kNN graph based on NN-Descent (with optional distance)
///
/// This function generates the kNN graph based via an approximate nearest
/// neighbour search based on the NN-Descent. The algorithm will use a
/// neighbours of neighbours logic to identify the approximate nearest
/// neighbours.
///
/// ### Params
///
/// * `mat` - Matrix in which rows represent the samples and columns the
///   respective embeddings for that sample
/// * `dist_metric` - The distance metric to use. One of `"euclidean"` or
///   `"cosine"`.
/// * `no_neighbours` - Number of neighbours for the KNN graph.
/// * `max_iter` - Maximum iterations for the algorithm.
/// * `delta` - Early stop criterium for the algorithm.
/// * `rho` - Sampling rate for the old neighbours. Will adaptively decrease
///   over time.
/// * `seed` - Seed for the NN Descent algorithm
/// * `verbose` - Controls verbosity of the algorithm
/// * `return_distances` - Shall the distances be returned.
///
/// ### Returns
///
/// The k-nearest neighbours based on the NN Desccent algorithm
///
/// ### Implementation details
///
/// In case of contrived synthetic data the algorithm sometimes does not
/// return enough neighbours. If that happens, the neighbours and distances will
/// be just padded.
#[allow(clippy::too_many_arguments)]
pub fn generate_knn_nndescent_with_dist(
    mat: MatRef<f32>,
    dist_metric: &str,
    no_neighbours: usize,
    max_iter: usize,
    delta: f32,
    rho: f32,
    seed: usize,
    verbose: bool,
    return_distances: bool,
) -> (Vec<Vec<usize>>, Option<Vec<Vec<f32>>>) {
    let graph: Vec<Vec<(usize, f32)>> = NNDescent::build(
        mat,
        no_neighbours,
        dist_metric,
        max_iter,
        delta,
        rho,
        seed,
        verbose,
    );

    let mut indices = Vec::with_capacity(graph.len());
    let mut distances = if return_distances {
        Some(Vec::with_capacity(graph.len()))
    } else {
        None
    };

    for (i, neighbours) in graph.into_iter().enumerate() {
        let mut ids: Vec<usize> = Vec::with_capacity(no_neighbours);
        let mut dists: Vec<f32> = Vec::with_capacity(no_neighbours);

        for (pid, dist) in neighbours {
            ids.push(pid);
            dists.push(dist);
        }

        if ids.len() < no_neighbours {
            let padding_needed = no_neighbours - ids.len();
            if ids.is_empty() {
                ids.resize(no_neighbours, i);
                dists.resize(no_neighbours, 0.0);
            } else {
                for j in 0..padding_needed {
                    ids.push(ids[j % ids.len()]);
                    dists.push(dists[j % dists.len()]);
                }
            }
        }

        indices.push(ids);
        if let Some(ref mut d) = distances {
            d.push(dists);
        }
    }

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
/// * `seed` - Seed for the Annoy algorithm
///
/// ### Returns
///
/// The k-nearest neighbours based on the Annoy algorithm
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

/// Get the kNN graph based on NN-Descent
///
/// This function generates the kNN graph based via an approximate nearest
/// neighbour search based on the NN-Descent. The algorithm will use a
/// neighbours of neighbours logic to identify the approximate nearest
/// neighbours.
///
/// ### Params
///
/// * `mat` - Matrix in which rows represent the samples and columns the
///   respective embeddings for that sample
/// * `dist_metric` - The distance metric to use. One of `"euclidean"` or
///   `"cosine"`.
/// * `no_neighbours` - Number of neighbours for the KNN graph.
/// * `max_iter` - Maximum iterations for the algorithm.
/// * `delta` - Early stop criterium for the algorithm.
/// * `rho` - Sampling rate for the old neighbours. Will adaptively decrease
///   over time.
/// * `seed` - Seed for the NN Descent algorithm
/// * `verbose` - Controls verbosity of the algorithm
///
/// ### Returns
///
/// The k-nearest neighbours based on the NN Desccent algorithm
///
/// ### Implementation details
///
/// In case of contrived synthetic data the algorithm sometimes does not
/// return enough neighbours. If that happens, the neighbours will be just
/// padded.
#[allow(clippy::too_many_arguments)]
pub fn generate_knn_nndescent(
    mat: MatRef<f32>,
    dist_metric: &str,
    no_neighbours: usize,
    max_iter: usize,
    delta: f32,
    rho: f32,
    seed: usize,
    verbose: bool,
) -> Vec<Vec<usize>> {
    let graph: Vec<Vec<(usize, f32)>> = NNDescent::build(
        mat,
        no_neighbours,
        dist_metric,
        max_iter,
        delta,
        rho,
        seed,
        verbose,
    );

    graph
        .into_iter()
        .enumerate()
        .map(|(i, neighbours)| {
            let mut ids: Vec<usize> = neighbours.into_iter().map(|(pid, _)| pid).collect();

            // pad if we don't have enough neighbours
            // need this to deal with the failing tests on weird synthetic data
            if ids.len() < no_neighbours {
                let padding_needed = no_neighbours - ids.len();
                if ids.is_empty() {
                    ids.resize(no_neighbours, i);
                } else {
                    for j in 0..padding_needed {
                        ids.push(ids[j % ids.len()]);
                    }
                }
            }

            ids
        })
        .collect()
}

/// Generate the kNN indices and distances
///
/// Helper function to generate kNN indices and distances in one go
///
/// ### Params
///
/// * `embd` - The embedding matrix to use to approximate neighbours and
///   calculate distances. Cells x features.
/// * `knn_params` - The parameters for the approximate nearest neighbour
///   search.
/// * `return_dist` - Return the distances.
/// * `seed` - Seed for reproducibility
/// * `verbose` - Controls verbosity of the function.
///
/// ### Returns
///
/// Tuple of `(indices of nearest neighbours, distances to these neighbours)`
pub fn generate_knn_with_dist(
    embd: MatRef<f32>,
    knn_params: &KnnParams,
    return_dist: bool,
    seed: usize,
    verbose: bool,
) -> (Vec<Vec<usize>>, Option<Vec<Vec<f32>>>) {
    // default to Annoy if nothing else...
    let knn_method = parse_knn_method(&knn_params.knn_method).unwrap_or(KnnSearch::Annoy);

    match knn_method {
        KnnSearch::Annoy => {
            let index_start = Instant::now();
            let annoy_index = build_annoy_index(embd, knn_params.n_tree, seed);
            let end_index = index_start.elapsed();

            if verbose {
                println!("Generated Annoy index: {:.2?}", end_index);
            }

            let query_start = Instant::now();
            let (mut indices, distances) = query_annoy_index(
                embd,
                &annoy_index,
                &knn_params.ann_dist,
                knn_params.k + 1, // Query k+1 to include self!
                knn_params.search_budget,
                return_dist,
                verbose,
            );
            let query_end = query_start.elapsed();

            if verbose {
                println!("Queried Annoy index: {:.2?}", query_end);
            }

            // remove first element (self) from each vector
            for idx_vec in indices.iter_mut() {
                idx_vec.remove(0);
            }

            let distances = distances.map(|mut dists| {
                for dist_vec in dists.iter_mut() {
                    dist_vec.remove(0);
                }
                dists
            });

            (indices, distances)
        }
        KnnSearch::Hnsw => {
            let index_start = Instant::now();
            let hnsw_index = build_hnsw_index(embd, &knn_params.ann_dist, seed);
            let end_index = index_start.elapsed();

            if verbose {
                println!("Generated HNSW index: {:.2?}", end_index);
            }

            let query_start = Instant::now();
            let (mut indices, distances) = query_hnsw_index(
                embd,
                &hnsw_index,
                &knn_params.ann_dist,
                knn_params.k + 1, // Query k+1 to include self!
                return_dist,
                verbose,
            );
            let query_end = query_start.elapsed();

            if verbose {
                println!("Queried HNSW index: {:.2?}", query_end);
            }

            // Remove first element (self) from each vector
            for idx_vec in indices.iter_mut() {
                idx_vec.remove(0);
            }

            let distances = distances.map(|mut dists| {
                for dist_vec in dists.iter_mut() {
                    dist_vec.remove(0);
                }
                dists
            });

            (indices, distances)
        }
        KnnSearch::NNDescent => generate_knn_nndescent_with_dist(
            embd,
            &knn_params.ann_dist,
            knn_params.k,
            knn_params.max_iter,
            knn_params.delta,
            knn_params.rho,
            seed,
            verbose,
            return_dist,
        ),
    }
}

///////////////////
// sNN functions //
///////////////////

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
