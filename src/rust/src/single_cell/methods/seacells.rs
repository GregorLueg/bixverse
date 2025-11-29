use extendr_api::List;
use faer::MatRef;
use rand::prelude::*;
use rand::rngs::StdRng;
use rayon::prelude::*;
use rustc_hash::FxHashSet;
use std::time::Instant;
use thousands::Separable;

use crate::core::base::sparse_math::*;
use crate::core::data::sparse_structures::*;
use crate::core::graph::knn::*;
use crate::single_cell::sc_knn_snn::KnnParams;

////////////
// Params //
////////////

/// Structure to store the SEACells parameters
///
/// ### Fields
///
/// **SEACells:**
///
/// * `n_sea_cells` - Number of sea cells to detect
/// * `max_fw_iters` - Maximum iterations for the Franke-Wolfe algorithm per
///   matrix update.
/// * `convergence_epsilon` - Defines the convergence threshold. Algorithm stops
///   when `RSS change < epsilon * RSS(0)`
/// * `max_iter` - Maximum iterations to run SEACells for
/// * `min_iter` - Minimum iterations to run SEACells for
/// * `prune_threshold` - The threshold below which values are set to 0 to
///   maintain sparsity and reduce memory pressure.
/// * `greedy_threshold` - Maximum number of cells, before defaulting to a more
///   rapid random selection of archetypes initially
/// * `pruning` - Shall tiny values during the Franke Wolfe updates be pruned.
///   This can affect numerical stability, but makes runs on large data sets
///   feasible.
/// * `pruning_threshold` - Values that should be pruned away.
///
/// **General kNN params**
///
/// * `k` - Number of neighbours for the kNN algorithm.
/// * `knn_method` - Which method to use for the generation of the kNN graph.
///   One of `"hnsw"`, `"annoy"` or `"nndescent"`
/// * `ann_dist` - The distance metric for the approximate nearest neighbour
///   search. For this method fixed to `"euclidean"`.
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
#[derive(Clone, Debug)]
pub struct SEACellsParams {
    // sea cell
    pub n_sea_cells: usize,
    pub max_fw_iters: usize,
    pub convergence_epsilon: f32,
    pub max_iter: usize,
    pub min_iter: usize,
    pub greedy_threshold: usize,
    pub graph_building: String,
    pub pruning: bool,
    pub pruning_threshold: f32,
    // general knn params
    pub knn_params: KnnParams,
}

impl SEACellsParams {
    /// Generate SEACellsParams from an R list
    ///
    /// Should values not be found within the List, the parameters will default
    /// to sensible defaults based on the original SEACells implementation.
    ///
    /// ### Params
    ///
    /// * `r_list` - The list with the SEACells parameters.
    ///
    /// ### Returns
    ///
    /// The `SEACellsParams` with all parameters set.
    pub fn from_r_list(r_list: List) -> Self {
        let knn_params = KnnParams::from_r_list(r_list.clone());

        let seacells_list = r_list.into_hashmap();

        let n_sea_cells = seacells_list
            .get("n_sea_cells")
            .and_then(|v| v.as_integer())
            .unwrap_or(0) as usize;

        let max_fw_iters = seacells_list
            .get("max_fw_iters")
            .and_then(|v| v.as_integer())
            .unwrap_or(50) as usize;

        // convergence_epsilon: algorithm converges when RSS change < epsilon * RSS(0)
        // Default: 1e-3 from Python implementation
        let convergence_epsilon = seacells_list
            .get("convergence_epsilon")
            .and_then(|v| v.as_real())
            .unwrap_or(1e-3) as f32;

        let max_iter = seacells_list
            .get("max_iter")
            .and_then(|v| v.as_integer())
            .unwrap_or(100) as usize;

        let min_iter = seacells_list
            .get("min_iter")
            .and_then(|v| v.as_integer())
            .unwrap_or(10) as usize;

        let greedy_threshold = seacells_list
            .get("greedy_threshold")
            .and_then(|v| v.as_integer())
            .unwrap_or(20000) as usize;

        let graph_building = seacells_list
            .get("graph_building")
            .and_then(|v| v.as_str())
            .unwrap_or("union")
            .to_string();

        let pruning = seacells_list
            .get("pruning")
            .and_then(|v| v.as_bool())
            .unwrap_or(false);

        let pruning_threshold: f32 = seacells_list
            .get("pruning_threshold")
            .and_then(|v| v.as_real())
            .unwrap_or(1e-7) as f32;

        Self {
            // seacell
            n_sea_cells,
            max_fw_iters,
            convergence_epsilon,
            max_iter,
            min_iter,
            greedy_threshold,
            graph_building,
            pruning,
            pruning_threshold,
            // knn
            knn_params,
        }
    }
}

/////////////
// Helpers //
/////////////

/// Convert SEACells hard assignments to metacell format
///
/// Transforms flat assignment vector (cell -> SEACell) into grouped format
/// (SEACell -> [cells]) suitable for aggregation functions.
///
/// ### Params
///
/// * `assignments` - Vector where assignments[cell_id] = seacell_id
/// * `k` - Number of SEACells
///
/// ### Returns
///
/// Vector of vectors, where result[seacell_id] contains all cells assigned to that SEACell
pub fn assignments_to_metacells(assignments: &[usize], k: usize) -> Vec<Vec<usize>> {
    let mut metacells = vec![Vec::new(); k];

    for (cell_id, &seacell_id) in assignments.iter().enumerate() {
        metacells[seacell_id].push(cell_id);
    }

    metacells
}

/// Convert sparse to dense with scaling in one pass
///
/// ### Params
///
/// * `mat` - The matrix to scale
/// * `scale` - The scale value
/// * `dense` - The slice to update
fn sparse_to_dense_csr_scaled(mat: &CompressedSparseData<f32>, scale: f32, dense: &mut [f32]) {
    let (nrows, ncols) = mat.shape;
    dense.fill(0.0);

    for row in 0..nrows {
        let row_start = mat.indptr[row];
        let row_end = mat.indptr[row + 1];

        for idx in row_start..row_end {
            let col = mat.indices[idx];
            dense[row * ncols + col] = mat.data[idx] * scale;
        }
    }
}

/// Helper function to prune tiny values and renormalise with L1
///
/// ### Params
///
/// * `mat` - Mutable reference to the CompressedSparseData to be pruned
/// * `threshold` - Pruning threshold
///
/// ### Returns
///
/// Pruned matrix.
fn prune_and_renormalise(mat: &mut CompressedSparseData<f32>, threshold: f32) {
    // Remove values below threshold
    let mut new_data = Vec::new();
    let mut new_indices = Vec::new();
    let mut new_indptr = vec![0];

    for row in 0..mat.shape.0 {
        let start = mat.indptr[row];
        let end = mat.indptr[row + 1];

        for idx in start..end {
            if mat.data[idx].abs() > threshold {
                new_data.push(mat.data[idx]);
                new_indices.push(mat.indices[idx]);
            }
        }
        new_indptr.push(new_data.len());
    }

    mat.data = new_data;
    mat.indices = new_indices;
    mat.indptr = new_indptr;

    // Renormalise columns to maintain sum-to-1 constraint
    normalise_csr_columns_l1(mat);
}

fn matrix_trace(mat: &CompressedSparseData<f32>) -> f32 {
    let n = mat.shape.0.min(mat.shape.1);
    let mut trace = 0.0;

    for i in 0..n {
        let row_start = mat.indptr[i];
        let row_end = mat.indptr[i + 1];

        for idx in row_start..row_end {
            if mat.indices[idx] == i {
                trace += mat.data[idx];
                break;
            }
        }
    }

    trace
}

/// Compute adaptive anisotropic diffusion kernel
///
/// Implementation from palantir package.  Uses knn/3-th nearest neighbor
/// distance as adaptive bandwidth. For edge (i,j) with distance d:
/// weight = exp(-d/σᵢ)
///
/// ### Params
///
/// * `knn_indices` - kNN indices for each cell
/// * `knn_distances` - kNN distances for each cell  
/// * `knn` - Number of nearest neighbours used
///
/// ### Returns
///
/// Symmetric kernel matrix
fn compute_diffusion_kernel(
    knn_indices: &[Vec<usize>],
    knn_distances: &[Vec<f32>],
    knn: usize,
) -> CompressedSparseData<f32> {
    let n = knn_indices.len();
    let adaptive_k = (knn / 3).max(1);

    let adaptive_std: Vec<f32> = knn_distances
        .iter()
        .map(|dists| {
            let mut sorted = dists.clone();
            sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
            sorted[adaptive_k - 1]
        })
        .collect();

    let mut rows = Vec::new();
    let mut cols = Vec::new();
    let mut vals = Vec::new();

    for (i, neighbours) in knn_indices.iter().enumerate() {
        for (idx, &j) in neighbours.iter().enumerate() {
            // need to square root here, as I am not doing this during kNN generation
            let dist = knn_distances[i][idx].sqrt();
            let weight = (-dist / adaptive_std[i]).exp();
            rows.push(i);
            cols.push(j);
            vals.push(weight);
        }
    }

    let w = coo_to_csr(&rows, &cols, &vals, (n, n));

    // symmetrise: kernel = W + W^T
    let w_t = w.transpose_and_convert();

    sparse_add_csr(&w, &w_t)
}

/// Compute diffusion maps from kernel matrix
///
/// Normalises kernel to transition matrix and performs eigendecomposition.
///
/// ### Params
///
/// * `kernel` - Symmetric kernel matrix
/// * `n_components` - Number of eigenvectors to compute
///
/// ### Returns
///
/// (eigenvalues, eigenvectors) where eigenvectors is (n × n_components)
fn diffusion_map_from_kernel(
    kernel: &mut CompressedSparseData<f32>,
    n_components: usize,
    seed: u64,
) -> (Vec<f32>, Vec<Vec<f32>>) {
    // Compute row sums (degrees)
    let row_sums: Vec<f32> = (0..kernel.shape.0)
        .map(|i| {
            (kernel.indptr[i]..kernel.indptr[i + 1])
                .map(|idx| kernel.data[idx])
                .sum()
        })
        .collect();

    // Symmetric normalisation: D^(-1/2) * K * D^(-1/2)
    for i in 0..kernel.shape.0 {
        let d_i_sqrt = row_sums[i].sqrt();
        for idx in kernel.indptr[i]..kernel.indptr[i + 1] {
            let j = kernel.indices[idx];
            let d_j_sqrt = row_sums[j].sqrt();
            kernel.data[idx] /= d_i_sqrt * d_j_sqrt;
        }
    }

    compute_largest_eigenpairs_lanczos(kernel, n_components, seed)
}

/// Determine multiscale space by scaling eigenvectors
///
/// Scales eigenvectors by λᵢ/(1-λᵢ) for diffusion distance metric.
///
/// ### Params
///
/// * `eigenvalues` - Eigenvalues from diffusion maps
/// * `eigenvectors` - Eigenvectors (n × n_components)
/// * `n_eigs` - Optional number of eigenvectors to use (None = auto-detect via
///   eigengap)
///
/// ### Returns
///
/// Scaled eigenvectors (n × n_eigs)
fn determine_multiscale_space(
    eigenvalues: &[f32],
    eigenvectors: &[Vec<f32>],
    n_eigs: Option<usize>,
) -> Vec<Vec<f32>> {
    let n = eigenvectors.len();

    // auto-detect n_eigs using eigengap if not provided
    let use_n_eigs = if let Some(n) = n_eigs {
        n
    } else {
        let gaps: Vec<f32> = eigenvalues.windows(2).map(|w| w[0] - w[1]).collect();

        let max_gap_idx = gaps
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .map(|(idx, _)| idx + 1)
            .unwrap_or(3);

        // DEBUG
        println!("Gaps: {:?}", &gaps[..gaps.len().min(10)]);
        println!("Max gap at index: {}", max_gap_idx);

        max_gap_idx.max(3).min(eigenvalues.len())
    };

    let use_indices: Vec<usize> = (1..use_n_eigs).collect();

    let mut scaled = vec![vec![0.0f32; use_indices.len()]; n];

    for (out_idx, &eig_idx) in use_indices.iter().enumerate() {
        let lambda = eigenvalues[eig_idx];
        let scale = lambda / (1.0 - lambda);

        for i in 0..n {
            scaled[i][out_idx] = eigenvectors[i][eig_idx] * scale;
        }
    }

    scaled
}

/// Max-min waypoint sampling
///
/// For each dimension, iteratively selects points maximizing the minimum
/// distance to already selected points.
///
/// ### Params
///
/// * `data` - Multiscale space (n × n_dims)
/// * `num_waypoints` - Target number of waypoints
/// * `seed` - Random seed for initial point selection
///
/// ### Returns
///
/// Indices of selected waypoints
fn max_min_sampling(data: &[Vec<f32>], num_waypoints: usize, seed: u64) -> Vec<usize> {
    let n = data.len();
    let n_dims = data[0].len();
    let no_iterations = (num_waypoints / n_dims).max(1);

    let mut rng = StdRng::seed_from_u64(seed);
    let mut waypoint_set = FxHashSet::default();

    for dim in 0..n_dims {
        let vec: Vec<f32> = data.iter().map(|row| row[dim]).collect();
        let mut iter_set = vec![rng.random_range(0..n)];
        let mut min_dists = vec![f32::MAX; n];

        // initialize distances to first point
        for i in 0..n {
            min_dists[i] = (vec[i] - vec[iter_set[0]]).abs();
        }

        // iteratively select maximally distant points
        for _ in 1..no_iterations {
            let new_wp = min_dists
                .iter()
                .enumerate()
                .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
                .map(|(idx, _)| idx)
                .unwrap();

            iter_set.push(new_wp);

            for i in 0..n {
                let dist = (vec[i] - vec[new_wp]).abs();
                min_dists[i] = min_dists[i].min(dist);
            }
        }

        waypoint_set.extend(iter_set);
    }

    waypoint_set.into_iter().collect()
}

//////////
// Main //
//////////

/// CPU implementation of the SEACells algorithm
///
/// SEACells identifies metacells (groupings of similar cells) using kernel
/// archetypal analysis. The algorithm solves a convex optimisation problem to
/// find archetypes that minimise reconstruction error whilst maintaining
/// sparsity.
///
/// This Rust implementation includes memory optimisations:
///
/// - Never materialises K_square to avoid O(n square) memory usage
/// - Prunes small values to maintain sparsity
/// - Supports fast random initialisation for large datasets
///
/// ### Fields
///
/// * `n_cells` - Number of cells in the dataset.
/// * `k` - Number of SEACells (metacells) to identify.
/// * `kernel_mat` - Sparse kernel matrix K computed from k-NN graph with RBF
///   weights.
/// * `a` - Assignment matrix (k × n) mapping cells to SEACells.
/// * `b` - Archetype matrix (n × k) defining SEACells as cell combinations.
/// * `archetypes` - Indices of cells selected as initial archetypes.
/// * `rss_history` - Residual sum of squares at each iteration.
/// * `convergence_threshold` - Absolute RSS change threshold for convergence.
/// * `params` - SEACell parameters.
pub struct SEACells<'a> {
    n_cells: usize,
    kernel_mat: Option<CompressedSparseData<f32>>,
    k_square: Option<CompressedSparseData<f32>>,
    a: Option<CompressedSparseData<f32>>,
    b: Option<CompressedSparseData<f32>>,
    archetypes: Option<Vec<usize>>,
    rss_history: Vec<f32>,
    convergence_threshold: Option<f32>,
    k_frobenius_norm_sq: Option<f32>,
    params: &'a SEACellsParams,
}

impl<'a> SEACells<'a> {
    /// Create a new SEACells instance
    ///
    /// ### Params
    ///
    /// * `n_cells` - Number of cells in the dataset
    /// * `n_seacells` - Number of SEACells to identify
    /// * `params` - Algorithm parameters
    ///
    /// ### Returns
    ///
    /// New `SEACells` instance with uninitialised matrices
    pub fn new(n_cells: usize, params: &'a SEACellsParams) -> Self {
        Self {
            n_cells,
            kernel_mat: None,
            k_square: None,
            a: None,
            b: None,
            archetypes: None,
            convergence_threshold: None,
            k_frobenius_norm_sq: None,
            rss_history: Vec::new(),
            params,
        }
    }

    /// Construct the kernel matrix from k-NN graph with adaptive RBF weights
    ///
    /// Builds a sparse kernel matrix K where K[i,j] represents similarity
    /// between cells i and j. Uses adaptive bandwidth RBF:
    ///
    /// ```exp(-dist^2/(σᵢ × σⱼ))```
    ///
    /// where σᵢ is the median distance to k nearest neighbours of cell i.
    ///
    /// The graph can be symmetrised using union (add edge if either direction
    /// exists) or intersection (add edge only if both directions exist).
    ///
    /// ### Params
    ///
    /// * `pca` - PCA/SVD matrix (n_cells × n_components)
    /// * `knn_indices` - k-NN indices for each cell
    /// * `knn_distances` - k-NN distances for each cell
    /// * `graph_construction` - Union or Intersection symmetrisation method
    /// * `verbose` - Print progress messages
    pub fn construct_kernel_mat(
        &mut self,
        pca: MatRef<f32>,
        knn_indices: &[Vec<usize>],
        knn_distances: &[Vec<f32>],
        verbose: bool,
    ) {
        let n = pca.nrows();
        let k = knn_indices[0].len();

        if verbose {
            println!("Computing adaptive bandwidth RBF kernel...");
        }

        let graph_construction = parse_knn_symmetrisation(&self.params.graph_building)
            .unwrap_or(KnnSymmetrisation::Intersection);

        let median_idx = k / 2;
        let median_dist = knn_distances
            .iter()
            .map(|d| d[median_idx].sqrt())
            .collect::<Vec<f32>>();

        let mut edges = FxHashSet::default();
        for (i, neighbours) in knn_indices.iter().enumerate() {
            for &j in neighbours {
                edges.insert((i, j));
            }
        }

        match graph_construction {
            KnnSymmetrisation::Union => {
                let to_add: Vec<_> = edges
                    .iter()
                    .filter_map(|&(i, j)| (!edges.contains(&(j, i))).then_some((j, i)))
                    .collect();
                edges.extend(to_add);
            }
            KnnSymmetrisation::Intersection => {
                let to_keep: FxHashSet<_> = edges
                    .iter()
                    .copied()
                    .filter(|&(i, j)| edges.contains(&(j, i)))
                    .collect();
                edges = to_keep;
            }
        }

        for i in 0..n {
            edges.insert((i, i));
        }

        let mut rows: Vec<usize> = Vec::new();
        let mut cols: Vec<usize> = Vec::new();
        let mut vals: Vec<f32> = Vec::new();

        for &(i, j) in &edges {
            let mut dist_square = 0_f32;
            for dim in 0..pca.ncols() {
                let diff = pca.get(i, dim) - pca.get(j, dim);
                dist_square += diff * diff;
            }
            let sigma_prod = median_dist[i] * median_dist[j];
            let val = (-dist_square / sigma_prod).exp();

            rows.push(i);
            cols.push(j);
            vals.push(val);
        }

        if verbose {
            println!(
                "Built kernel with {} non-zeros",
                vals.len().separate_with_underscores()
            );
        }

        let kernel = coo_to_csr(&rows, &cols, &vals, (n, n));
        let k_square = csr_matmul_csr_optimised(&kernel, &kernel.transpose_and_convert());

        if verbose {
            println!(
                "K_square has {} non-zeros",
                k_square.data.len().separate_with_underscores()
            );
        }

        if self.n_cells > 20000 {
            if verbose {
                println!("Pre-computing kernel Frobenius norm due to large data set...");
            }
            let k_frob = frobenius_norm(&kernel);
            self.k_frobenius_norm_sq = Some(k_frob * k_frob);
        }

        self.kernel_mat = Some(kernel);
        self.k_square = Some(k_square);
    }

    /// Fit the SEACells model
    ///
    /// Runs the main optimisation loop:
    ///
    /// 1. Initialises archetypes (greedy CSSP or random)
    /// 2. Initialises A and B matrices
    /// 3. Alternates updating A and B using Frank-Wolfe until convergence
    ///
    /// Convergence is reached when RSS change < epsilon × RSS(0), subject to
    /// minimum iteration requirements.
    ///
    /// ### Params
    ///
    /// * `seed` - Random seed for reproducibility
    /// * `verbose` - Print progress and RSS values
    pub fn fit(&mut self, seed: usize, verbose: bool) {
        assert!(
            self.kernel_mat.is_some(),
            "Must construct kernel matrix first"
        );
        assert!(self.archetypes.is_some(), "Must find archetypes first");

        self.initialise_matrices(verbose, seed as u64);

        let a = self.a.as_ref().unwrap();
        let b = self.b.as_ref().unwrap();

        let initial_rss = self.compute_rss(a, b);
        self.rss_history.push(initial_rss);
        self.convergence_threshold = Some(self.params.convergence_epsilon * initial_rss);

        if verbose {
            println!("Initial RSS: {:.6}", initial_rss);
            println!(
                "Convergence threshold: {:.6}",
                self.convergence_threshold.unwrap()
            );
        }

        let mut converged = false;
        let mut n_iter = 0;

        while (!converged && n_iter < self.params.max_iter) || n_iter < self.params.min_iter {
            let iter_start = Instant::now();
            n_iter += 1;

            let b_current = self.b.clone().unwrap();
            let a_current = self.a.clone().unwrap();

            let a_new = self.update_a_mat(&b_current, &a_current, verbose);
            let b_new = self.update_b_mat(&a_new, &b_current, verbose);

            let rss = self.compute_rss(&a_new, &b_new);
            self.rss_history.push(rss);

            self.a = Some(a_new);
            self.b = Some(b_new);

            let iter_duration = iter_start.elapsed();

            if verbose {
                println!(
                    "Iteration {}: RSS = {:.6}, Time = {:.2}s",
                    n_iter,
                    rss,
                    iter_duration.as_secs_f32()
                );
            }

            if n_iter > 1 {
                let rss_diff = (self.rss_history[n_iter - 1] - self.rss_history[n_iter]).abs();
                if rss_diff < self.convergence_threshold.unwrap() {
                    if n_iter >= self.params.min_iter {
                        if verbose {
                            println!("Converged after {} iterations!", n_iter);
                        }
                        converged = true;
                    } else if verbose && !converged {
                        println!(
                            "Convergence criteria met at iteration {} (continuing to min_iter)",
                            n_iter
                        );
                        converged = true;
                    }
                }
            }
        }

        if !converged && verbose {
            println!(
                "Warning: Algorithm did not converge after {} iterations",
                self.params.max_iter
            );
        }
    }

    /// Initialise archetypes using adaptive strategy
    ///
    /// For small datasets (< greedy_threshold): combines waypoint + greedy CSSP
    /// For large datasets (>= greedy_threshold): uses fast random initialisation
    ///
    /// ### Params
    ///
    /// * `knn_indices` - k-NN indices for each cell
    /// * `knn_distances` - k-NN distances for each cell
    /// * `verbose` - Print which method is selected
    /// * `seed` - Random seed for initialisation
    pub fn initialise_archetypes(
        &mut self,
        knn_indices: &[Vec<usize>],
        knn_distances: &[Vec<f32>],
        verbose: bool,
        seed: u64,
    ) {
        if self.n_cells > self.params.greedy_threshold {
            if verbose {
                println!(
                    "Dataset large (n={}), using fast random init (threshold: {})",
                    self.n_cells.separate_with_underscores(),
                    self.params.greedy_threshold
                );
            }
            self.initialise_archetypes_random(verbose, seed);
        } else {
            self.initialise_archetypes_combined(knn_indices, knn_distances, verbose, seed);
        }
    }

    /// Fast random archetype initialisation
    ///
    /// Randomly samples k cells as initial archetypes. Used for large datasets
    /// where greedy CSSP is computationally expensive.
    ///
    /// ### Params
    ///
    /// * `verbose` - Print number of archetypes selected
    /// * `seed` - Random seed for reproducibility
    fn initialise_archetypes_random(&mut self, verbose: bool, seed: u64) {
        let mut rng = StdRng::seed_from_u64(seed);
        let mut indices: Vec<usize> = (0..self.n_cells).collect();
        indices.shuffle(&mut rng);

        let archetypes: Vec<usize> = indices.into_iter().take(self.params.n_sea_cells).collect();

        if verbose {
            println!("Selected {} random archetypes", archetypes.len());
        }

        self.archetypes = Some(archetypes);
    }

    /// Combined waypoint + greedy initialisation (matches Python logic)
    ///
    /// 1. Gets waypoint centres (may return < k cells)
    /// 2. Tops up with greedy CSSP to reach k cells
    /// 3. Deduplicates and takes first k unique cells
    ///
    /// ### Params
    ///
    /// * `knn_indices` - k-NN indices for each cell
    /// * `knn_distances` - k-NN distances for each cell
    /// * `verbose` - Print selection counts
    /// * `seed` - Random seed for waypoint sampling
    fn initialise_archetypes_combined(
        &mut self,
        knn_indices: &[Vec<usize>],
        knn_distances: &[Vec<f32>],
        verbose: bool,
        seed: u64,
    ) {
        let k = self.params.n_sea_cells;

        // get waypoint centres
        if verbose {
            println!("Computing diffusion maps for waypoint initialisation...");
        }

        let mut kernel =
            compute_diffusion_kernel(knn_indices, knn_distances, self.params.knn_params.k);

        let (eigenvalues, eigenvectors) =
            diffusion_map_from_kernel(&mut kernel, self.params.knn_params.k, seed);

        let multiscale = determine_multiscale_space(&eigenvalues, &eigenvectors, Some(10));

        let waypoint_ix = max_min_sampling(&multiscale, k, seed);

        if verbose {
            println!(
                "Selecting {} cells from waypoint initialisation.",
                waypoint_ix.len()
            );
        }

        // calculate how many more needed from greedy
        let from_greedy = k.saturating_sub(waypoint_ix.len());

        // get greedy centres with +10 buffer (for deduplication)
        if verbose {
            println!("Initialising residual matrix using greedy column selection");
        }
        let greedy_ix = self.get_greedy_centres(from_greedy + 10);

        if verbose {
            println!(
                "Selecting {} cells from greedy initialisation.",
                from_greedy
            );
        }

        // combine and deduplicate
        let mut all_ix = waypoint_ix;
        all_ix.extend(greedy_ix);

        let mut seen = FxHashSet::default();
        let unique_ix: Vec<usize> = all_ix
            .into_iter()
            .filter(|&x| seen.insert(x))
            .take(k)
            .collect();

        self.archetypes = Some(unique_ix);
    }

    /// Get greedy centres (extracted from initialise_archetypes_greedy)
    ///
    /// ### Params
    ///
    /// * `n_centres` - Number of centres to select
    ///
    /// ### Returns
    ///
    /// Vector of selected cell indices
    fn get_greedy_centres(&self, n_centres: usize) -> Vec<usize> {
        let k_square = self.k_square.as_ref().unwrap();
        let n = k_square.shape.0;

        let results: Vec<(Vec<f32>, f32)> = (0..n)
            .into_par_iter()
            .map(|i| {
                let mut e_i = vec![0_f32; n];
                e_i[i] = 1.0;
                let k2_e_i = csr_matvec(k_square, &e_i);
                let g_val = k2_e_i[i];
                (k2_e_i, g_val)
            })
            .collect();

        let mut f = vec![0_f32; n];
        let mut g = vec![0_f32; n];
        for (i, (k2_e_i, g_val)) in results.iter().enumerate() {
            for j in 0..n {
                f[j] += k2_e_i[j] * k2_e_i[j];
            }
            g[i] = *g_val;
        }

        let mut omega: Vec<Vec<f32>> = vec![vec![0_f32; n]; n_centres];
        let mut centres: Vec<usize> = Vec::with_capacity(n_centres);

        let mut e_p = vec![0.0f32; n];
        let mut omega_new = vec![0.0f32; n];

        for iter in 0..n_centres {
            let mut best_idx = 0;
            let mut best_score = f32::MIN;

            for i in 0..n {
                if g[i] > 1e-15 {
                    let score = f[i] / g[i];
                    if score > best_score {
                        best_score = score;
                        best_idx = i;
                    }
                }
            }

            centres.push(best_idx);

            e_p.fill(0.0);
            e_p[best_idx] = 1.0;
            let k2_col = csr_matvec(k_square, &e_p);

            let mut delta = k2_col.clone();
            for i in 0..n {
                let omega_sum: f32 = (0..iter).map(|r| omega[r][best_idx] * omega[r][i]).sum();
                delta[i] -= omega_sum;
            }

            delta[best_idx] = delta[best_idx].max(0.0);
            let delta_p_sqrt = delta[best_idx].sqrt().max(1e-6);

            for i in 0..n {
                omega_new[i] = delta[i] / delta_p_sqrt;
            }

            let omega_sq_norm: f32 = omega_new.iter().map(|&x| x * x).sum();
            let k_omega_new = csr_matvec(k_square, &omega_new);

            for i in 0..n {
                let omega_hadamard = omega_new[i] * omega_new[i];
                let term1 = omega_sq_norm * omega_hadamard;

                let pl: f32 = (0..iter).map(|r| omega[r][best_idx] * omega[r][i]).sum();
                let term2 = omega_new[i] * (k_omega_new[i] - pl);

                f[i] += -2.0 * term2 + term1;
                g[i] += omega_hadamard;
            }

            omega[iter].copy_from_slice(&omega_new);
        }

        centres
    }

    /// Initialise A and B matrices
    ///
    /// Creates:
    ///
    /// - B matrix: one-hot encoding of archetype cells (n × k)
    /// - A matrix: random sparse assignments normalised to sum to 1 per cell
    ///   (k × n)
    ///
    /// Each cell is randomly assigned to ~25% of archetypes with random weights,
    /// then normalised. A is then updated once using Frank-Wolfe for better
    /// starting point.
    ///
    /// ### Params
    ///
    /// * `verbose` - Print initialisation message
    /// * `seed` - Random seed for A matrix initialisation
    fn initialise_matrices(&mut self, verbose: bool, seed: u64) {
        let archetypes = self.archetypes.as_ref().unwrap();
        let k = archetypes.len();
        let n = self.n_cells;

        if verbose {
            println!("Initialising A and B matrices...");
        }

        let mut b_rows = Vec::new();
        let mut b_cols = Vec::new();
        let mut b_vals = Vec::new();

        for (col, &row) in archetypes.iter().enumerate() {
            b_rows.push(row);
            b_cols.push(col);
            b_vals.push(1_f32);
        }

        let b = coo_to_csr(&b_rows, &b_cols, &b_vals, (n, k));

        let archetypes_per_cell = (k as f32 * 0.25).ceil() as usize;
        let mut rng = StdRng::seed_from_u64(seed);

        let mut a_rows = Vec::new();
        let mut a_cols = Vec::new();
        let mut a_vals = Vec::new();

        for cell in 0..n {
            for _ in 0..archetypes_per_cell {
                let archetype = rng.random_range(0..k);
                a_rows.push(archetype);
                a_cols.push(cell);
                a_vals.push(rng.random::<f32>());
            }
        }

        let mut a = coo_to_csr(&a_rows, &a_cols, &a_vals, (k, n));

        normalise_csr_columns_l1(&mut a);

        a = self.update_a_mat(&b, &a, verbose);

        self.a = Some(a);
        self.b = Some(b);
    }

    /// Update assignment matrix A using Frank-Wolfe algorithm
    ///
    /// Solves:
    ///
    /// ```min ||K_square - K_square @ B @ A||^2```
    ///
    /// subject to A columns summing to 1.
    ///
    /// Computes gradient G = 2(t1 @ A - t2) where:
    /// - t1 = (K_square @ B)^T @ B
    /// - t2 = (K_square @ B)^T
    ///
    /// For each cell, sets weight to 1 for archetype with minimum gradient,
    /// then takes convex step towards this solution.
    ///
    /// ### Params
    ///
    /// * `b` - Current archetype matrix
    /// * `a_prev` - Previous assignment matrix
    ///
    /// ### Returns
    ///
    /// Updated assignment matrix
    fn update_a_mat(
        &self,
        b: &CompressedSparseData<f32>,
        a_prev: &CompressedSparseData<f32>,
        verbose: bool,
    ) -> CompressedSparseData<f32> {
        let k_square = self.k_square.as_ref().unwrap();

        let k_b = csr_matmul_csr_optimised(k_square, b);
        let t2 = k_b.transpose_and_convert();
        let t1 = csr_matmul_csr_optimised(&t2, b);

        let mut a = a_prev.clone();
        let n = a.shape.1;
        let k = a.shape.0;

        let mut g_dense = vec![0.0f32; k * n];
        let mut e_data = Vec::with_capacity(n);

        for t in 0..self.params.max_fw_iters {
            let t1_a = csr_matmul_csr_optimised(&t1, &a);
            let g_mat = sparse_subtract_csr(&t1_a, &t2);

            sparse_to_dense_csr_scaled(&g_mat, 2.0, &mut g_dense);

            e_data.clear();
            let argmins: Vec<usize> = (0..n)
                .into_par_iter()
                .map(|col| {
                    let mut min_val = g_dense[col];
                    let mut min_idx = 0;

                    for row in 1..k {
                        let val = g_dense[row * n + col];
                        if val < min_val {
                            min_val = val;
                            min_idx = row;
                        }
                    }
                    min_idx
                })
                .collect();

            for (col, &row) in argmins.iter().enumerate() {
                e_data.push((row, col, 1.0f32));
            }
            e_data.sort_unstable_by_key(|&(r, c, _)| (r, c));

            let e_rows: Vec<usize> = e_data.iter().map(|&(r, _, _)| r).collect();
            let e_cols: Vec<usize> = e_data.iter().map(|&(_, c, _)| c).collect();
            let e_vals: Vec<f32> = e_data.iter().map(|&(_, _, v)| v).collect();

            let e = coo_to_csr_presorted(&e_rows, &e_cols, &e_vals, (k, n));

            let step_size = 2.0 / (t as f32 + 2.0);

            let e_minus_a = sparse_subtract_csr(&e, &a);
            let update = sparse_scalar_multiply_csr(&e_minus_a, step_size);
            a = sparse_add_csr(&a, &update);

            if self.params.pruning {
                prune_and_renormalise(&mut a, self.params.pruning_threshold);
            }

            if verbose && (t + 1) % 10 == 0 {
                println!(
                    "  A matrix Frank-Wolfe iteration: {} / {}",
                    t + 1,
                    self.params.max_fw_iters
                );
            }
        }

        a
    }

    /// Update archetype matrix B using Frank-Wolfe algorithm
    ///
    /// Solves:
    ///
    /// ```min ||K_square - K_square @ B @ A||^2```
    ///
    /// subject to B columns summing to 1.
    ///
    /// Computes gradient G = 2(K_square @ B @ t1 - t2) where:
    ///
    /// - t1 = A @ A^T
    /// - t2 = K_square @ A^T
    ///
    /// For each archetype, sets weight to 1 for cell with minimum gradient,
    /// then takes convex step towards this solution.
    ///
    /// ### Params
    ///
    /// * `a` - Current assignment matrix
    /// * `b_prev` - Previous archetype matrix
    ///
    /// ### Returns
    ///
    /// Updated archetype matrix
    fn update_b_mat(
        &self,
        a: &CompressedSparseData<f32>,
        b_prev: &CompressedSparseData<f32>,
        verbose: bool,
    ) -> CompressedSparseData<f32> {
        let k_square = self.k_square.as_ref().unwrap();

        let a_t = a.transpose_and_convert();
        let t1 = csr_matmul_csr_optimised(a, &a_t);
        let t2 = csr_matmul_csr_optimised(k_square, &a_t);

        // Early stopping threshold
        const FW_TOLERANCE: f32 = 1e-4;

        let mut b = b_prev.clone();
        let n = b.shape.0;
        let k = b.shape.1;

        let mut g_dense = vec![0.0f32; n * k];
        let mut e_data = Vec::with_capacity(k);

        for t in 0..self.params.max_fw_iters {
            let k_b = csr_matmul_csr_optimised(k_square, &b);
            let k_b_t1 = csr_matmul_csr_optimised(&k_b, &t1);
            let g_mat = sparse_subtract_csr(&k_b_t1, &t2);

            sparse_to_dense_csr_scaled(&g_mat, 2.0, &mut g_dense);

            e_data.clear();
            let argmins: Vec<usize> = (0..k)
                .into_par_iter()
                .map(|col| {
                    let mut min_val = g_dense[col];
                    let mut min_idx = 0;

                    for row in 1..n {
                        let val = g_dense[row * k + col];
                        if val < min_val {
                            min_val = val;
                            min_idx = row;
                        }
                    }
                    min_idx
                })
                .collect();

            for (col, &row) in argmins.iter().enumerate() {
                e_data.push((row, col, 1.0f32));
            }
            e_data.sort_unstable_by_key(|&(r, c, _)| (r, c));

            let e_rows: Vec<usize> = e_data.iter().map(|&(r, _, _)| r).collect();
            let e_cols: Vec<usize> = e_data.iter().map(|&(_, c, _)| c).collect();
            let e_vals: Vec<f32> = e_data.iter().map(|&(_, _, v)| v).collect();

            let e = coo_to_csr_presorted(&e_rows, &e_cols, &e_vals, (n, k));

            let step_size = 2.0 / (t as f32 + 2.0);

            let e_minus_b = sparse_subtract_csr(&e, &b);
            let update = sparse_scalar_multiply_csr(&e_minus_b, step_size);

            // for early stops
            let update_norm: f32 = update.data.iter().map(|&x| x * x).sum::<f32>().sqrt();

            b = sparse_add_csr(&b, &update);

            if self.params.pruning {
                prune_and_renormalise(&mut b, self.params.pruning_threshold);
            }

            if verbose && (t + 1) % 10 == 0 {
                println!(
                    "  B matrix Frank-Wolfe iteration: {} / {}",
                    t + 1,
                    self.params.max_fw_iters
                );
            }

            // early stopping if update is negligible
            if update_norm < FW_TOLERANCE && t >= 10 {
                if verbose {
                    println!("  B matrix FW converged early at iteration {}", t + 1);
                }
                break;
            }
        }

        b
    }

    /// Compute residual sum of squares (RSS)
    ///
    /// Calculates Frobenius norm: ```||K - K @ B @ A||_F^2```
    ///
    /// Note: Uses K (not K²) for reconstruction to measure approximation quality.
    ///
    /// ### Params
    ///
    /// * `a` - Assignment matrix
    /// * `b` - Archetype matrix
    ///
    /// ### Returns
    ///
    /// RSS value (lower is better fit)
    fn compute_rss(&self, a: &CompressedSparseData<f32>, b: &CompressedSparseData<f32>) -> f32 {
        // Use simple method for small datasets, trace method for large ones
        if self.n_cells <= 20000 {
            self.compute_rss_simple(a, b)
        } else {
            self.compute_rss_trace(a, b)
        }
    }

    /// Fast RSS computation for small datasets (materialises reconstruction)
    ///
    /// This version is quite fast and works well on small data sets.
    ///
    /// ### Params
    ///
    /// * `a` - The A matrix
    /// * `b` - The B matrix
    ///
    /// ### Returns
    ///
    /// The residual sum of squares (RSS)
    fn compute_rss_simple(
        &self,
        a: &CompressedSparseData<f32>,
        b: &CompressedSparseData<f32>,
    ) -> f32 {
        let k_mat = self.kernel_mat.as_ref().unwrap();
        let k_b = csr_matmul_csr_optimised(k_mat, b);
        let reconstruction = csr_matmul_csr_optimised(&k_b, a);
        let diff = sparse_subtract_csr(k_mat, &reconstruction);
        frobenius_norm(&diff)
    }

    /// Memory-efficient RSS computation for large datasets (uses trace trick)
    ///
    /// This version is slower, but does not blow up memory.
    ///
    /// ### Params
    ///
    /// * `a` - The A matrix
    /// * `b` - The B matrix
    ///
    /// ### Returns
    ///
    /// The residual sum of squares (RSS)
    fn compute_rss_trace(
        &self,
        a: &CompressedSparseData<f32>,
        b: &CompressedSparseData<f32>,
    ) -> f32 {
        let k_square = self.k_square.as_ref().unwrap();

        // Term 1: ||K||_F^2 (cached)
        let k_frob_sq = self.k_frobenius_norm_sq.unwrap();

        // Term 2: -2 * trace(K_square @ B @ A)
        // Reorder: trace(K_square @ B @ A) = trace(A @ K_square @ B)
        let k2_b = csr_matmul_csr_optimised(k_square, b); // (n×n) @ (n×k) = (n×k)
        let a_k2b = csr_matmul_csr_optimised(a, &k2_b); // (k×n) @ (n×k) = (k×k) -> SMALL!
        let trace_term = matrix_trace(&a_k2b);

        // Term 3: trace(A^T @ B^T @ K_square @ B @ A)
        // Reorder: trace(A @ A^T @ B^T @ K_square @ B)
        let a_t = a.transpose_and_convert(); // (n×k)
        let a_at = csr_matmul_csr_optimised(a, &a_t); // (k×n) @ (n×k) = (k×k) -> SMALL!

        let b_t = b.transpose_and_convert(); // (k×n)
        let bt_k2b = csr_matmul_csr_optimised(&b_t, &k2_b); // (k×n) @ (n×k) = (k×k) -> SMALL!

        let result = csr_matmul_csr_optimised(&a_at, &bt_k2b); // (k×k) @ (k×k) = (k×k) -> SMALL!
        let reconstruction_frob_sq = matrix_trace(&result);

        (k_frob_sq - 2.0 * trace_term + reconstruction_frob_sq).sqrt()
    }

    /// Get hard cell assignments (each cell assigned to one SEACell)
    ///
    /// ### Returns
    ///
    /// Vector of SEACell assignments (0 to k-1)
    pub fn get_hard_assignments(&self) -> Vec<usize> {
        let a = self.a.as_ref().expect("Model not fitted yet");
        let n = a.shape.1;
        let k = a.shape.0;

        let mut assignments = vec![0usize; n];

        for cell in 0..n {
            let mut max_val = f32::NEG_INFINITY;
            let mut max_idx = 0;

            for archetype in 0..k {
                let row_start = a.indptr[archetype];
                let row_end = a.indptr[archetype + 1];

                for idx in row_start..row_end {
                    if a.indices[idx] == cell {
                        if a.data[idx] > max_val {
                            max_val = a.data[idx];
                            max_idx = archetype;
                        }
                        break;
                    }
                }
            }

            assignments[cell] = max_idx;
        }

        assignments
    }

    /// Get RSS history
    ///
    /// ### Returns
    ///
    /// Vectors of the RSS
    pub fn get_rss_history(&self) -> &[f32] {
        &self.rss_history
    }

    /// Get the archetype cell indices
    ///
    /// ### Returns
    ///
    /// Vector of cell indices selected as archetypes
    ///
    /// ### Panics
    ///
    /// Panics if archetypes have not been initialised yet
    pub fn get_archetypes(&self) -> Vec<usize> {
        self.archetypes
            .as_ref()
            .expect("Archetypes not initialised yet")
            .clone()
    }
}
