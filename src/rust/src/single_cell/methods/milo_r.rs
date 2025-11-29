use extendr_api::List;
use faer::MatRef;
use instant_distance::{HnswMap, Search};
use rayon::prelude::*;

use crate::core::base::stats::median;
use crate::core::graph::annoy::*;
use crate::core::graph::knn::*;
use crate::single_cell::sc_knn_snn::*;

/////////////
// Helpers //
/////////////

/// Structure for MiloR algorithm parameters
///
/// ### Fields
///
/// **MiloR params**
///
/// * `prop` - Proportion of cells to sample as neighbourhood indices
/// * `k_refine` - Number of neighbours to use for refinement
/// * `refinement_strategy` - Strategy for refining sampled indices
///   (`"approximate"`, `"bruteforce"`, or `"index"`)
/// * `index_type` - Type of kNN index to use (`"annoy"` or `"hnsw"`)
///
/// **General kNN params**
///
/// * `knn_params` - All of the kNN parameters
pub struct MiloRParams {
    // MiloR parameters
    pub prop: f64,
    pub k_refine: usize,
    pub refinement_strategy: String,
    pub index_type: String,
    // kNN
    pub knn_params: KnnParams,
}

impl MiloRParams {
    /// Generate MiloRParams from an R list
    ///
    /// Should values not be found within the List, the parameters will default
    /// to sensible defaults based on heuristics.
    ///
    /// ### Params
    ///
    /// * `r_list` - The list with the MiloR parameters
    ///
    /// ### Returns
    ///
    /// The `MiloRParams` with all parameters set.
    pub fn from_r_list(r_list: List) -> Self {
        let knn_params = KnnParams::from_r_list(r_list.clone());

        let params_list = r_list.into_hashmap();

        let prop = params_list
            .get("prop")
            .and_then(|v| v.as_real())
            .unwrap_or(0.2);

        let k_refine = params_list
            .get("k_refine")
            .and_then(|v| v.as_integer())
            .unwrap_or(20) as usize;

        let index_type = std::string::String::from(
            params_list
                .get("index_type")
                .and_then(|v| v.as_str())
                .unwrap_or("annoy"),
        );

        let refinement_strategy = std::string::String::from(
            params_list
                .get("refinement_strategy")
                .and_then(|v| v.as_str())
                .unwrap_or("approximate"),
        );

        Self {
            prop,
            k_refine,
            refinement_strategy,
            index_type,
            knn_params,
        }
    }
}

/// Enum wrapper for different kNN index implementations
///
/// ### Variants
///
/// * `Annoy` - Approximate nearest neighbour index using trees
/// * `Hnsw` - Hierarchical navigable small world graph index
pub enum KnnIndex {
    /// The Annoy index
    Annoy(AnnoyIndex),
    /// The HNSW index
    Hnsw(HnswMap<Point, usize>),
}

impl KnnIndex {
    /// Generate a new instance of the kNN index
    ///
    /// ### Params
    ///
    /// * `embd` - The embedding matrix of cells x features to use to the
    ///   the generation
    /// * `knn_params` - The KnnParams with distance type, number of trees, etc.
    /// * `seed` - Random seed for reproducibility
    ///
    /// ### Returns
    ///
    /// Initialised `KnnIndex`.
    pub fn new(
        embd: MatRef<f32>,
        index_type: KnnIndexType,
        knn_params: &KnnParams,
        seed: usize,
    ) -> Self {
        match index_type {
            KnnIndexType::AnnoyIndex => {
                KnnIndex::Annoy(AnnoyIndex::new(embd, knn_params.n_tree, seed))
            }
            KnnIndexType::HnswIndex => {
                KnnIndex::Hnsw(build_hnsw_index(embd, &knn_params.ann_dist, seed))
            }
        }
    }

    /// Query for k nearest neighbours of a single point
    ///
    /// ### Params
    ///
    /// * `query_point` - The slice of values defining the query point
    /// * `knn_params` - The KnnParams with distance type, search budget, etc.
    /// * `k` - Number of neighbours to return
    ///
    /// ### Returns
    ///
    /// Tuple of `(neighbour indices, distances to neighbours)`
    pub fn query_single(
        &self,
        query_point: &[f32],
        knn_params: &KnnParams,
        k: usize,
    ) -> (Vec<usize>, Vec<f32>) {
        match self {
            KnnIndex::Annoy(index) => {
                let metric = parse_ann_dist(&knn_params.ann_dist).unwrap_or(AnnDist::Euclidean);

                index.query(query_point, &metric, k, Some(knn_params.search_budget))
            }
            KnnIndex::Hnsw(index) => {
                let metric = parse_ann_dist(&knn_params.ann_dist).unwrap_or(AnnDist::Euclidean);
                let point = Point(query_point.to_vec(), metric);
                let mut search = Search::default();

                let results: Vec<_> = index.search(&point, &mut search).take(k).collect();

                let indices = results.iter().map(|item| *item.value).collect();
                let distances = results.iter().map(|item| item.distance).collect();

                (indices, distances)
            }
        }
    }
}

/// Enum specifying which kNN index type to use
///
/// ### Variants
///
/// * `AnnoyIndex` - Use Annoy index
/// * `HnswIndex` - Use HNSW index
pub enum KnnIndexType {
    /// Annoy
    AnnoyIndex,
    /// HNSW
    HnswIndex,
}

//////////////
// Sampling //
//////////////

/// Enum specifying the refinement strategy for neighbourhood sampling
///
/// ### Variants
///
/// * `Approximate` - Search within k neighbours only
/// * `BruteForce` - Linear search through all cells
/// * `IndexBased` - Use existing kNN index for search
#[derive(Debug, Clone, Copy)]
pub enum RefinementStrategy {
    Approximate,
    BruteForce,
    IndexBased,
}

/// Helper function to parse the refinement strategy
///
/// ### Params
///
/// * `s` - String specifying the strategy to use
///
/// ### Returns
///
/// The Option of the chosen `RefinementStrategy`
pub fn parse_refinement_strategy(s: &str) -> Option<RefinementStrategy> {
    match s.to_lowercase().as_str() {
        "approximate" => Some(RefinementStrategy::Approximate),
        "bruteforce" => Some(RefinementStrategy::BruteForce),
        "index" => Some(RefinementStrategy::IndexBased),
        _ => None,
    }
}

/// Helper function to parse the kNN index type
///
/// ### Params
///
/// * `s` - String specifying which kNN index to use
///
/// ### Returns
///
/// The Option of the chosen `KnnIndexType`
pub fn parse_index_type(s: &str) -> Option<KnnIndexType> {
    match s.to_lowercase().as_str() {
        "annoy" => Some(KnnIndexType::AnnoyIndex),
        "hnsw" => Some(KnnIndexType::HnswIndex),
        _ => None,
    }
}

/// Helper function to compute the median positions
///
/// ### Params
///
/// * `embd` - The embedding matrix that was used for the generation of the kNN
///   graph.
/// * `neighbours` - Slice of indices for the neighbours
///
/// ### Returns
///
/// Vector of median features
fn compute_median_position(embd: MatRef<f32>, neighbours: &[usize]) -> Vec<f32> {
    let n_feature = embd.ncols();

    let mut median_point = vec![0.0f32; n_feature];
    for feat_idx in 0..n_feature {
        let values = neighbours
            .iter()
            .map(|&nb_idx| embd[(nb_idx, feat_idx)])
            .collect::<Vec<f32>>();

        median_point[feat_idx] = median(&values).unwrap_or(0_f32);
    }

    median_point
}

/// Find the cell nearest to a median position within a subset of candidates
///
/// ### Params
///
/// * `embd` - The embedding matrix
/// * `median_point` - The median position to query
/// * `candidates` - Indices of candidate cells to search within
/// * `metric` - Distance metric to use
///
/// ### Returns
///
/// Index of the nearest cell within the candidate subset
fn find_nearest_in_subset(
    embd: MatRef<f32>,
    median_point: &[f32],
    candidates: &[usize],
    metric: &AnnDist,
) -> usize {
    let median_row = MatRef::from_row_major_slice(median_point, 1, embd.ncols()).row(0);

    candidates
        .par_iter()
        .map(|&idx| {
            let dist = compute_distance_knn(median_row, embd.row(idx), metric);
            (idx, dist)
        })
        .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
        .unwrap()
        .0
}

/// Find the cell nearest to a median position using brute force search
///
/// ### Params
///
/// * `embd` - The embedding matrix
/// * `median_point` - The median position to query
/// * `metric` - Distance metric to use
///
/// ### Returns
///
/// Index of the nearest cell in the entire dataset
fn find_nearest_bruteforce(embd: MatRef<f32>, median_point: &[f32], metric: &AnnDist) -> usize {
    let median_row = MatRef::from_row_major_slice(median_point, 1, embd.ncols()).row(0);

    (0..embd.nrows())
        .into_par_iter()
        .map(|idx| {
            let dist = compute_distance_knn(median_row, embd.row(idx), metric);
            (idx, dist)
        })
        .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
        .unwrap()
        .0
}

/// Find the cell nearest to a median position using the kNN index
///
/// ### Params
///
/// * `index` - The kNN index to query
/// * `knn_params` - Parameters for the kNN search
/// * `median_point` - The median position to query
///
/// ### Returns
///
/// Index of the nearest cell
fn find_nearest_with_index(
    index: &KnnIndex,
    knn_params: &KnnParams,
    median_point: &[f32],
) -> usize {
    let (indices, _) = index.query_single(median_point, knn_params, 1);

    indices[0]
}

/// Refine neighbourhood sampling by shifting indices towards local median
/// positions
///
/// ### Params
///
/// * `embd` - The embedding matrix
/// * `knn_indices` - The kNN graph as adjacency list
/// * `sampled_indices` - Initial sampled cell indices
/// * `k_refine` - Number of neighbours to use for computing median
/// * `knn_params` - Parameters for distance calculation
/// * `strategy` - Refinement strategy to use
/// * `knn_index` - Optional kNN index for index-based strategy
/// * `verbose` - Whether to print progress messages
///
/// ### Returns
///
/// Refined cell indices after shifting to nearest median positions
#[allow(clippy::too_many_arguments)]
pub fn refine_sampling_with_strategy(
    embd: MatRef<f32>,
    knn_indices: &[Vec<usize>],
    sampled_indices: &[usize],
    k_refine: usize,
    knn_params: &KnnParams,
    strategy: &RefinementStrategy,
    knn_index: Option<&KnnIndex>,
    verbose: bool,
) -> Vec<usize> {
    if verbose {
        println!("Running refined sampling");
    }

    let mut refined = Vec::with_capacity(sampled_indices.len());

    let dist_metric = parse_ann_dist(&knn_params.ann_dist).unwrap_or(AnnDist::Euclidean);

    for &sample_idx in sampled_indices {
        let mut neighbours = Vec::with_capacity(k_refine);
        for j in 0..k_refine.min(knn_indices[0].len()) {
            let neighbour_idx = knn_indices[sample_idx][j];
            neighbours.push(neighbour_idx);
        }

        let median_point = compute_median_position(embd, &neighbours);

        let best_idx = match strategy {
            RefinementStrategy::Approximate => {
                find_nearest_in_subset(embd, &median_point, &neighbours, &dist_metric)
            }
            RefinementStrategy::BruteForce => {
                find_nearest_bruteforce(embd, &median_point, &dist_metric)
            }
            RefinementStrategy::IndexBased => {
                if let Some(index) = knn_index {
                    find_nearest_with_index(index, knn_params, &median_point)
                } else {
                    // Fallback to brute force
                    find_nearest_bruteforce(embd, &median_point, &dist_metric)
                }
            }
        };

        refined.push(best_idx);
    }

    refined
}

/// Compute distances to the k-th nearest neighbour for each index cell
///
/// ### Params
///
/// * `embd` - The embedding matrix
/// * `knn_indices` - The kNN graph as adjacency list
/// * `index_cells` - Indices of neighbourhood centre cells
/// * `kth_col` - Which neighbour to compute distance to (0-indexed)
///
/// ### Returns
///
/// Vector of distances to k-th neighbour for each index cell
pub fn compute_kth_distances_from_matrix(
    embd: MatRef<f32>,
    knn_indices: &[Vec<usize>],
    index_cells: &[usize],
    kth_col: usize,
) -> Vec<f64> {
    index_cells
        .par_iter()
        .map(|&cell_idx| {
            let kth_neighbour = knn_indices[cell_idx][kth_col];

            compute_distance_knn(
                embd.row(cell_idx),
                embd.row(kth_neighbour),
                &AnnDist::Euclidean,
            ) as f64
        })
        .collect()
}

/// Build sparse neighbourhood matrix in COO (triplet) format
///
/// Each neighbourhood includes the index cell plus its k nearest neighbours.
///
/// ### Params
///
/// * `knn_indices` - The kNN graph as adjacency list
/// * `index_cells` - Indices of neighbourhood centre cells
///
/// ### Returns
///
/// Tuple of `(row_indices, col_indices, values)` in COO format where
/// each neighbourhood is a column and non-zero entries indicate membership
pub fn build_nhood_matrix(
    knn_indices: &[Vec<usize>],
    index_cells: &[usize],
) -> (Vec<usize>, Vec<usize>, Vec<f64>) {
    let k = knn_indices[0].len();
    let n_nhoods = index_cells.len();

    // Pre-allocate (over-estimate)
    let mut row_indices = Vec::with_capacity(n_nhoods * (k + 1));
    let mut col_indices = Vec::with_capacity(n_nhoods * (k + 1));
    let mut values = Vec::with_capacity(n_nhoods * (k + 1));

    for (nh_idx, &cell_idx) in index_cells.iter().enumerate() {
        row_indices.push(cell_idx);
        col_indices.push(nh_idx);
        values.push(1.0);

        for j in 0..k {
            let neighbor_idx = knn_indices[cell_idx][j];
            row_indices.push(neighbor_idx);
            col_indices.push(nh_idx);
            values.push(1.0);
        }
    }

    (row_indices, col_indices, values)
}
