use ann_search_rs::utils::KnnValidation;
use ann_search_rs::*;
use extendr_api::List;
use faer::{MatRef, RowRef};
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use std::time::Instant;

use crate::core::data::sparse_structures::*;
use crate::core::graph::knn::AnnDist;

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
#[derive(Default)]
pub enum KnnSearch {
    /// Annoy-based
    #[default]
    Annoy,
    /// Hierarchical Navigable Small World
    Hnsw,
    /// NNDescent
    NNDescent,
    /// Exhaustive
    Exhaustive,
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
/// * `k` - Number of neighbours to search.
///
/// **Annoy**
///
/// * `n_tree` - Number of trees for the generation of the index.
/// * `search_budget` - Optional search budget. If not provided, will default
///   to `k * n_trees * 20`. Good ranges for the multipler are 2 to 20.
///
/// **NN Descent**
///
/// * `delta` - Early termination criterium.
/// * `diversify_prob` - Diversifying probability at the end of the index
///   generation.
/// * `ef_budget` - Optional query budget.
///
/// **LSH**
///
/// * `bits` - Number of bits to use.
/// * `n_tables` - Number of hash tables to use.
/// * `max_candidates` - Optional query budget.
///
/// **IVF**
///
/// * `n_centroids` - Number of centroids to use.
/// * `n_probes` - Number of centroids to probe.
#[derive(Clone, Debug)]
pub struct KnnParams {
    // general params
    pub knn_method: String,
    pub ann_dist: String,
    pub k: usize,
    // annoy params
    pub n_tree: usize,
    pub search_budget: Option<usize>,
    // nn descent params
    pub diversify_prob: f32,
    pub delta: f32,
    pub ef_budget: Option<usize>,
    // hnsw
    pub m: usize,
    pub ef_construction: usize,
    pub ef_search: usize,
}

impl KnnParams {
    /// Generate KnnParams from an R list
    ///
    /// Should values not be found within the List, the parameters will default
    /// to sensible defaults based on heuristics.
    ///
    /// ### Params
    ///
    /// * `r_list` - The list with the kNN parameters.
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
            .get("n_trees")
            .and_then(|v| v.as_integer())
            .unwrap_or(50) as usize;

        let search_budget = params_list
            .get("search_budget")
            .and_then(|v| v.as_integer())
            .map(|v| v as usize);

        // nn descent
        let diversify_prob = params_list
            .get("diversify_prob")
            .and_then(|v| v.as_real())
            .unwrap_or(0.0) as f32;

        let delta = params_list
            .get("delta")
            .and_then(|v| v.as_real())
            .unwrap_or(0.001) as f32;

        let ef_budget = params_list
            .get("ef_budget")
            .and_then(|v| v.as_integer())
            .map(|v| v as usize);

        // hnsw
        let m = params_list
            .get("m")
            .and_then(|v| v.as_integer())
            .unwrap_or(16) as usize;

        let ef_construction = params_list
            .get("ef_construction")
            .and_then(|v| v.as_integer())
            .unwrap_or(200) as usize;

        let ef_search = params_list
            .get("ef_search")
            .and_then(|v| v.as_integer())
            .unwrap_or(100) as usize;

        Self {
            knn_method,
            ann_dist,
            k,
            n_tree,
            search_budget,
            ef_budget,
            diversify_prob,
            delta,
            m,
            ef_construction,
            ef_search,
        }
    }

    /// Generate a version of this with sensible base parameters
    ///
    /// ### Returns
    ///
    /// Self.
    pub fn new() -> Self {
        Self {
            knn_method: "annoy".to_string(),
            ann_dist: "cosine".to_string(),
            k: 15,
            n_tree: 50,
            search_budget: None,
            diversify_prob: 0.5,
            delta: 0.001,
            ef_budget: None,
            m: 16,
            ef_construction: 200,
            ef_search: 100,
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
        "exhaustive" => Some(KnnSearch::Exhaustive),
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

////////////////////
// Main functions //
////////////////////

/// Helper function to abstract out common patterns
///
/// ### Params
///
/// * `no_neighbours` - Number of neighbours
/// * `seed` - Seed for reproducibility
/// * `verbose` - Controls verbosity of the function
/// * `build_index` - Build index function
/// * `query_index` - Query index self
/// * `validate_index` - Self validation of the data
/// * `index_name` - Name of the index
///
/// ### Returns
///
/// The kNN graph
fn build_and_query_knn<I>(
    no_neighbours: usize,
    verbose: bool,
    build_index: impl FnOnce() -> I,
    query_index: impl FnOnce(&I) -> (Vec<Vec<usize>>, Option<Vec<Vec<f32>>>),
    index_name: &str,
) -> (Vec<Vec<usize>>, I) {
    let start = Instant::now();
    let index = build_index();
    if verbose {
        println!("Generated {} index: {:.2?}", index_name, start.elapsed());
    }

    let start = Instant::now();
    let (indices, _) = query_index(&index);

    let res: Vec<Vec<usize>> = indices
        .into_iter()
        .enumerate()
        .map(|(i, mut neighbors)| {
            neighbors.retain(|&x| x != i);
            neighbors.truncate(no_neighbours);
            neighbors
        })
        .collect();

    if verbose {
        println!(
            "Identified approximate nearest neighbours via {}: {:.2?}",
            index_name,
            start.elapsed()
        );
    }

    (res, index)
}

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
/// * `m` - Number of connections per layer (M parameter)
/// * `ef_const` - Size of dynamic candidate list during construction
/// * `ef_search` - Size of candidate list during search (higher = better
///   recall, slower)
/// * `seed` - Seed for the HNSW algorithm
/// * `verbose` - Controls verbosity
///
/// ### Returns
///
/// The k-nearest neighbours based on the HNSW algorithm. Function does not
/// return self.
#[allow(clippy::too_many_arguments)]
pub fn generate_knn_hnsw(
    mat: MatRef<f32>,
    dist_metric: &str,
    no_neighbours: usize,
    m: usize,
    ef_const: usize,
    ef_search: usize,
    seed: usize,
    verbose: bool,
) -> Vec<Vec<usize>> {
    let (res, index) = build_and_query_knn(
        no_neighbours,
        verbose,
        || build_hnsw_index(mat, m, ef_const, dist_metric, seed, verbose),
        |idx| query_hnsw_self(idx, no_neighbours + 1, ef_search, false, true),
        "HNSW",
    );

    if verbose {
        let recall = index.validate_index(no_neighbours, seed, None);
        println!(
            "Recall of approximate nearest neighbours search in random subset: {:.2}",
            recall
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
/// * `search_budget` - Optional search budget per given query. If not provided,
///   it will use `k * n_trees * 20`.
/// * `seed` - Seed for the Annoy algorithm
///
/// ### Returns
///
/// The k-nearest neighbours based on the Annoy algorithm. Function does not
/// return self.
pub fn generate_knn_annoy(
    mat: MatRef<f32>,
    dist_metric: &str,
    no_neighbours: usize,
    n_trees: usize,
    search_budget: Option<usize>,
    seed: usize,
    verbose: bool,
) -> Vec<Vec<usize>> {
    let (res, index) = build_and_query_knn(
        no_neighbours,
        verbose,
        || build_annoy_index(mat, dist_metric.to_string(), n_trees, seed),
        |idx| query_annoy_self(idx, no_neighbours + 1, search_budget, false, verbose),
        "Annoy",
    );

    if verbose {
        let recall = index.validate_index(no_neighbours, seed, None);
        println!(
            "Recall of approximate nearest neighbours search in random subset: {:.2}",
            recall
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
/// * `diversify_prob` - How many of the edges in the index shall be diversified
///   after index generation.
/// * `ef_budget` - Optional query search budget.
/// * `delta` - Early stop criterium for the algorithm.
/// * `seed` - Seed for the NN Descent algorithm
/// * `verbose` - Controls verbosity of the algorithm
///
/// ### Returns
///
/// The k-nearest neighbours based on the NNDescent algorithm. Function does not
/// return self.
#[allow(clippy::too_many_arguments)]
pub fn generate_knn_nndescent(
    mat: MatRef<f32>,
    dist_metric: &str,
    no_neighbours: usize,
    diversify_prob: f32,
    ef_budget: Option<usize>,
    delta: f32,
    seed: usize,
    verbose: bool,
) -> Vec<Vec<usize>> {
    let (res, index) = build_and_query_knn(
        no_neighbours,
        verbose,
        || {
            build_nndescent_index(
                mat,
                dist_metric,
                delta,
                diversify_prob,
                None,
                None,
                None,
                None,
                seed,
                verbose,
            )
        },
        |idx| query_nndescent_self(idx, no_neighbours + 1, ef_budget, false, verbose),
        "NNDescent",
    );

    if verbose {
        let recall = index.validate_index(no_neighbours, seed, None);
        println!(
            "Recall of approximate nearest neighbours search in random subset: {:.2}",
            recall
        );
    }

    res
}

/// Get the kNN graph based on an exhaustive search
///
/// ### Params
///
/// * `mat` - Matrix in which rows represent the samples and columns the
///   respective embeddings for that sample
/// * `dist_metric` - The distance metric to use. One of `"euclidean"` or
///   `"cosine"`.
/// * `no_neighbours` - Number of neighbours for the KNN graph.
/// * `verbose` - Controls verbosity of the function
///
/// ### Returns
///
/// The k-nearest neighbours based on the exhaustive linear search. Function
/// does not return self.
pub fn generate_knn_exhaustive(
    mat: MatRef<f32>,
    dist_metric: &str,
    no_neighbours: usize,
    verbose: bool,
) -> Vec<Vec<usize>> {
    let (res, _) = build_and_query_knn(
        no_neighbours,
        verbose,
        || build_exhaustive_index(mat, dist_metric),
        |idx| query_exhaustive_self(idx, no_neighbours + 1, false, verbose),
        "exhaustive linear search",
    );
    res
}

///////////////////
// With distance //
///////////////////

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
    // first helper function to remove self
    fn remove_self(
        mut indices: Vec<Vec<usize>>,
        distances: Option<Vec<Vec<f32>>>,
    ) -> (Vec<Vec<usize>>, Option<Vec<Vec<f32>>>) {
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

    // second helper function to time everything
    fn timed<T>(name: &str, verbose: bool, f: impl FnOnce() -> T) -> T {
        let start = Instant::now();
        let result = f();
        if verbose {
            println!("{}: {:.2?}", name, start.elapsed());
        }
        result
    }

    let knn_method = parse_knn_method(&knn_params.knn_method).unwrap_or_default();
    let k_plus_one = knn_params.k + 1;

    let (indices, distances) = match knn_method {
        KnnSearch::Annoy => {
            let index = timed("Generated Annoy index", verbose, || {
                build_annoy_index(embd, knn_params.ann_dist.clone(), knn_params.n_tree, seed)
            });
            let (indices, distances) = timed("Queried Annoy index", verbose, || {
                query_annoy_index(
                    embd,
                    &index,
                    k_plus_one,
                    knn_params.search_budget,
                    return_dist,
                    verbose,
                )
            });
            if verbose {
                let recall = index.validate_index(k_plus_one, seed, None);
                println!(
                    "Recall of approximate nearest neighbours search in random subset: {:.2}",
                    recall
                );
            }
            (indices, distances)
        }
        KnnSearch::Hnsw => {
            let index = timed("Generated HNSW index", verbose, || {
                build_hnsw_index(
                    embd,
                    knn_params.m,
                    knn_params.ef_construction,
                    &knn_params.ann_dist,
                    seed,
                    verbose,
                )
            });
            let (indices, distances) = timed("Queried HNSW index", verbose, || {
                query_hnsw_index(
                    embd,
                    &index,
                    k_plus_one,
                    knn_params.ef_search,
                    return_dist,
                    verbose,
                )
            });
            if verbose {
                let recall = index.validate_index(k_plus_one, seed, None);
                println!(
                    "Recall of approximate nearest neighbours search in random subset: {:.2}",
                    recall
                );
            }
            (indices, distances)
        }
        KnnSearch::NNDescent => {
            let index = timed("Generated NNDescent index", verbose, || {
                build_nndescent_index(
                    embd,
                    &knn_params.ann_dist,
                    knn_params.delta,
                    knn_params.diversify_prob,
                    None,
                    None,
                    None,
                    None,
                    seed,
                    verbose,
                )
            });
            let (indices, distances) = timed("Queried NNDescent index", verbose, || {
                query_nndescent_index(
                    embd,
                    &index,
                    k_plus_one,
                    knn_params.ef_budget,
                    true,
                    verbose,
                )
            });
            if verbose {
                let recall = index.validate_index(k_plus_one, seed, None);
                println!(
                    "Recall of approximate nearest neighbours search in random subset: {:.2}",
                    recall
                );
            }
            (indices, distances)
        }
        KnnSearch::Exhaustive => {
            let index = timed("Generated Exhaustive index", verbose, || {
                build_exhaustive_index(embd, &knn_params.knn_method)
            });
            timed("Queried Exhaustive index", verbose, || {
                query_exhaustive_index(embd, &index, k_plus_one, true, verbose)
            })
        }
    };

    remove_self(indices, distances)
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
