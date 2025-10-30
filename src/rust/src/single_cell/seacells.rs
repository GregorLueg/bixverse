use extendr_api::List;
use faer::MatRef;
use rand::prelude::*;
use rand::rngs::StdRng;
use rustc_hash::FxHashSet;
use std::time::Instant;

use crate::core::data::sparse_structures::*;

///////////
// Enums //
///////////

/// SNN similarity method
#[derive(Clone, Copy)]
pub enum SeaCellGraphGen {
    /// Only intersecting nearest neigbhbours will be considered
    Intersection,
    /// The union of nearest neighbours will be considered
    Union,
}

/// Helper function to parse the SEACell graph generation
///
/// ### Params
///
/// * `s` - Type of graph to build
///
/// ### Returns
///
/// Option of the SeaCellGraphGen
pub fn parse_seacell_graph(s: &str) -> Option<SeaCellGraphGen> {
    match s.to_lowercase().as_str() {
        "intersection" => Some(SeaCellGraphGen::Intersection),
        "union" => Some(SeaCellGraphGen::Union),
        _ => None,
    }
}

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
///
/// **General kNN params**
///
/// * `k` - Number of neighbours for the kNN algorithm.
/// * `knn_method` - Which method to use for the generation of the kNN graph.
///   One of `"hnsw"`, `"annoy"` or `"nndescent"`
/// * `ann_dist` - The distance metric for the approximate nearest neighbour
///   search. One of `"cosine"` or `"euclidean"`.
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
    // general knn params
    pub k: usize,
    pub knn_method: String,
    pub ann_dist: String,
    // annoy params
    pub n_trees: usize,
    pub search_budget: usize,
    // nn descent params
    pub nn_max_iter: usize,
    pub rho: f32,
    pub delta: f32,
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

        // generall knn
        let k = seacells_list
            .get("k")
            .and_then(|v| v.as_integer())
            .unwrap_or(25) as usize;
        let knn_method = seacells_list
            .get("knn_method")
            .and_then(|v| v.as_str())
            .unwrap_or("annoy")
            .to_string();
        let ann_dist = seacells_list
            .get("ann_dist")
            .and_then(|v| v.as_str())
            .unwrap_or("cosine")
            .to_string();

        // annoy
        let n_trees = seacells_list
            .get("n_trees")
            .and_then(|v| v.as_integer())
            .unwrap_or(100) as usize;

        let search_budget = seacells_list
            .get("search_budget")
            .and_then(|v| v.as_integer())
            .unwrap_or(100) as usize;

        // nn descent
        let nn_max_iter = seacells_list
            .get("nn_max_iter")
            .and_then(|v| v.as_integer())
            .unwrap_or(15) as usize;

        let rho = seacells_list
            .get("rho")
            .and_then(|v| v.as_real())
            .unwrap_or(1.0) as f32;

        let delta = seacells_list
            .get("delta")
            .and_then(|v| v.as_real())
            .unwrap_or(0.001) as f32;

        Self {
            // seacell
            n_sea_cells,
            max_fw_iters,
            convergence_epsilon,
            max_iter,
            min_iter,
            greedy_threshold,
            graph_building,
            // knn
            k,
            knn_method,
            ann_dist,
            n_trees,
            search_budget,
            nn_max_iter,
            rho,
            delta,
        }
    }
}

/////////////
// Helpers //
/////////////

/// Helper function to calculate the kernel matrix with a vector
///
/// This is WAY faster and memory-efficient than creating the full kernel
/// matrix in memory...
///
/// ### Params
///
/// * `kernel` - The kernel matrix
/// * `v` - The vector to multiply the matrix with
///
/// ### Returns
///
/// The final vector
fn kernel_squared_matvec(kernel: &CompressedSparseData<f32>, v: &[f32]) -> Vec<f32> {
    let kernel_t = kernel.transpose_and_convert();
    let temp = csr_matvec(&kernel_t, v);

    csr_matvec(kernel, &temp)
}

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

/// Convert CSR sparse matrix to dense row-major vector
///
/// ### Params
///
/// * `mat` - Sparse matrix in CSR format
///
/// ### Returns
///
/// Dense vector in row-major order (row_idx * ncols + col_idx)
fn sparse_to_dense_csr(mat: &CompressedSparseData<f32>) -> Vec<f32> {
    let (nrows, ncols) = mat.shape;
    let mut dense = vec![0.0f32; nrows * ncols];

    for row in 0..nrows {
        let row_start = mat.indptr[row];
        let row_end = mat.indptr[row + 1];

        for idx in row_start..row_end {
            let col = mat.indices[idx];
            dense[row * ncols + col] = mat.data[idx];
        }
    }

    dense
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

        let graph_construction = parse_seacell_graph(&self.params.graph_building)
            .unwrap_or(SeaCellGraphGen::Intersection);

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
            SeaCellGraphGen::Union => {
                let to_add: Vec<_> = edges
                    .iter()
                    .filter_map(|&(i, j)| (!edges.contains(&(j, i))).then_some((j, i)))
                    .collect();
                edges.extend(to_add);
            }
            SeaCellGraphGen::Intersection => {
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
            println!("Built kernel with {} non-zeros", vals.len());
        }

        let kernel = coo_to_csr(&rows, &cols, &vals, (n, n));
        let k_square = csr_matmul_csr(&kernel, &kernel.transpose_and_convert());

        if verbose {
            println!("K_square has {} non-zeros", k_square.data.len());
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

        self.initialise_archetypes(verbose, seed as u64);
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

            let a_new = self.update_a_mat(&b_current, &a_current);
            let b_new = self.update_b_mat(&a_new, &b_current);

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
                    if verbose {
                        println!("Converged after {} iterations!", n_iter);
                    }
                    converged = true;
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
    /// Selects between greedy CSSP (accurate but slow) and random selection
    /// (fast but less optimal) based on dataset size threshold.
    ///
    /// ### Params
    ///
    /// * `verbose` - Print which method is selected
    /// * `seed` - Random seed for random initialisation
    fn initialise_archetypes(&mut self, verbose: bool, seed: u64) {
        if self.n_cells > self.params.greedy_threshold {
            if verbose {
                println!(
                    "Dataset large (n={}), using fast random init (threshold: {})",
                    self.n_cells, self.params.greedy_threshold
                );
            }
            self.initialise_archetypes_random(verbose, seed);
        } else {
            if verbose {
                println!("Dataset small (n={}), using greedy CSSP", self.n_cells);
            }
            self.initialise_archetypes_greedy(verbose);
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

    /// Greedy CSSP archetype initialisation
    ///
    /// Uses Column Subset Selection Problem (CSSP) to select archetypes that
    /// maximally cover the data manifold. Iteratively selects cells that best
    /// explain remaining variance in K_square.
    ///
    /// Uses on-the-fly K_square computation to avoid materialising the full
    /// matrix. Complexity: O(k × n) K_square column computations.
    ///
    /// ### Params
    ///
    /// * `verbose` - Print progress every 10 archetypes
    fn initialise_archetypes_greedy(&mut self, verbose: bool) {
        let kernel = self.kernel_mat.as_ref().unwrap();
        let n = kernel.shape.0;
        let k = self.params.n_sea_cells;

        if verbose {
            println!("Initialising {} archetypes via greedy CSSP...", k);
        }

        let mut f = vec![0_f32; n];
        let mut g = vec![0_f32; n];

        for i in 0..n {
            let mut e_i = vec![0_f32; n];
            e_i[i] = 1.0;

            let k2_e_i = kernel_squared_matvec(kernel, &e_i);

            for j in 0..n {
                f[j] += k2_e_i[j] * k2_e_i[j];
            }
            g[i] = k2_e_i[i];
        }

        let mut omega: Vec<Vec<f32>> = vec![vec![0_f32; n]; k];
        let mut centers: Vec<usize> = Vec::with_capacity(k);

        for iter in 0..k {
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

            centers.push(best_idx);

            let mut e_p = vec![0.0f32; n];
            e_p[best_idx] = 1.0;
            let k2_col = kernel_squared_matvec(kernel, &e_p);

            let mut delta = k2_col.clone();
            for i in 0..n {
                let mut omega_sum = 0.0f32;
                for r in 0..iter {
                    omega_sum += omega[r][best_idx] * omega[r][i];
                }
                delta[i] -= omega_sum;
            }

            delta[best_idx] = delta[best_idx].max(0.0);

            let delta_p_sqrt = delta[best_idx].sqrt().max(1e-6);
            let mut omega_new = vec![0.0f32; n];
            for i in 0..n {
                omega_new[i] = delta[i] / delta_p_sqrt;
            }

            let omega_sq_norm: f32 = omega_new.iter().map(|&x| x * x).sum();

            let k_omega_new = kernel_squared_matvec(kernel, &omega_new);

            let dots: Vec<f32> = (0..iter)
                .map(|r| omega[r].iter().zip(&omega_new).map(|(a, b)| a * b).sum())
                .collect();

            for i in 0..n {
                let omega_hadamard = omega_new[i] * omega_new[i];
                let term1 = omega_sq_norm * omega_hadamard;

                let pl: f32 = (0..iter).map(|r| dots[r] * omega[r][i]).sum();
                let term2 = omega_new[i] * (k_omega_new[i] - pl);

                f[i] += -2.0 * term2 + term1;
                g[i] += omega_hadamard;
            }

            omega[iter] = omega_new;

            if verbose && (iter + 1) % 10 == 0 {
                println!("  Selected {} / {} archetypes", iter + 1, k);
            }
        }

        self.archetypes = Some(centers);
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

        a = self.update_a_mat(&b, &a);

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
    ) -> CompressedSparseData<f32> {
        let k_square = self.k_square.as_ref().unwrap();

        let k_b = csr_matmul_csr(k_square, b);
        let t2 = k_b.transpose_and_convert();
        let t1 = csr_matmul_csr(&t2, b);

        let mut a = a_prev.clone();
        let n = a.shape.1;
        let k = a.shape.0;

        for t in 0..self.params.max_fw_iters {
            let t1_a = csr_matmul_csr(&t1, &a);
            let g_mat = sparse_subtract_csr(&t1_a, &t2);
            let g_scaled = sparse_scalar_multiply_csr(&g_mat, 2.0);

            let g_dense = sparse_to_dense_csr(&g_scaled);

            let mut e_rows = Vec::with_capacity(n);
            let mut e_cols = Vec::with_capacity(n);
            let mut e_vals = Vec::with_capacity(n);

            for col in 0..n {
                let mut min_val = g_dense[col];
                let mut min_idx = 0;

                for row in 1..k {
                    let val = g_dense[row * n + col];
                    if val < min_val {
                        min_val = val;
                        min_idx = row;
                    }
                }

                e_rows.push(min_idx);
                e_cols.push(col);
                e_vals.push(1.0f32);
            }

            let e = coo_to_csr(&e_rows, &e_cols, &e_vals, (k, n));

            let step_size = 2.0 / (t as f32 + 2.0);
            let e_minus_a = sparse_subtract_csr(&e, &a);
            let update = sparse_scalar_multiply_csr(&e_minus_a, step_size);
            a = sparse_add_csr(&a, &update);
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
    ) -> CompressedSparseData<f32> {
        let k_square = self.k_square.as_ref().unwrap();

        let a_t = a.transpose_and_convert();
        let t1 = csr_matmul_csr(a, &a_t);
        let t2 = csr_matmul_csr(k_square, &a_t);

        let mut b = b_prev.clone();
        let n = b.shape.0;
        let k = b.shape.1;

        for t in 0..self.params.max_fw_iters {
            let k_b = csr_matmul_csr(k_square, &b);
            let k_b_t1 = csr_matmul_csr(&k_b, &t1);

            let g_mat = sparse_subtract_csr(&k_b_t1, &t2);
            let g_scaled = sparse_scalar_multiply_csr(&g_mat, 2.0);

            let g_dense = sparse_to_dense_csr(&g_scaled);

            let mut e_rows = Vec::with_capacity(k);
            let mut e_cols = Vec::with_capacity(k);
            let mut e_vals = Vec::with_capacity(k);

            for col in 0..k {
                let mut min_val = g_dense[col];
                let mut min_idx = 0;

                for row in 1..n {
                    let val = g_dense[row * k + col];
                    if val < min_val {
                        min_val = val;
                        min_idx = row;
                    }
                }

                e_rows.push(min_idx);
                e_cols.push(col);
                e_vals.push(1.0f32);
            }

            let e = coo_to_csr(&e_rows, &e_cols, &e_vals, (n, k));

            let step_size = 2.0 / (t as f32 + 2.0);
            let e_minus_b = sparse_subtract_csr(&e, &b);
            let update = sparse_scalar_multiply_csr(&e_minus_b, step_size);
            b = sparse_add_csr(&b, &update);
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
        let k_mat = self.kernel_mat.as_ref().unwrap();

        let k_b = csr_matmul_csr(k_mat, b);
        let reconstruction = csr_matmul_csr(&k_b, a);

        let diff = sparse_subtract_csr(k_mat, &reconstruction);

        frobenius_norm(&diff)
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
}
