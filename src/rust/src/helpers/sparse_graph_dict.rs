use extendr_api::prelude::*;

use faer::{
    linalg::solvers::{PartialPivLu, Solve},
    ColRef, Mat, MatRef,
};
use rand::{rngs::StdRng, Rng, SeedableRng};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::time::Instant;

use crate::helpers::graph::{adjacency_to_laplacian, get_knn_graph_adj};
use crate::helpers::linalg::{column_cosine, sylvester_solver};

////////////////////////
// Params and results //
////////////////////////

/// Structure to store the dual graph regularised dictionary learning params
///
/// ### Fields
///
/// * `sparsity` - Sparsity constraint (max non-zero coefficients per signal)
/// * `dict_size` - Dictionary size
/// * `alpha` - Sample context regularisation weight
/// * `beta` - Feature effect regularisation weight
/// * `max_iter` - Maximum number of iterations
/// * `k_neighbours` - Number of neighbours in the KNN graph
/// * `admm_iter` - ADMM iterations for sparse coding
/// * `rho` - ADMM step size
#[derive(Debug, Clone, Default)]
pub struct DgrdlParams {
    pub sparsity: usize,
    pub dict_size: usize,
    pub alpha: f64,
    pub beta: f64,
    pub max_iter: usize,
    pub k_neighbours: usize,
    pub admm_iter: usize,
    pub rho: f64,
}

#[allow(dead_code)]
impl DgrdlParams {
    /// Generate the DGRDL parameters from an R list
    ///
    /// If values are not found, will use default values
    ///
    /// ### Params
    ///
    /// * `r_list` - The R list containing the parameters
    ///
    /// ### Returns
    ///
    /// The `DgrdlParams` structure based on the R list
    pub fn from_r_list(r_list: List) -> Self {
        let dgrdl_params = r_list.into_hashmap();

        let sparsity = dgrdl_params
            .get("sparsity")
            .and_then(|v| v.as_integer())
            .unwrap_or(5) as usize;
        let dict_size = dgrdl_params
            .get("dict_size")
            .and_then(|v| v.as_integer())
            .unwrap_or(5) as usize;
        let alpha = dgrdl_params
            .get("alpha")
            .and_then(|v| v.as_real())
            .unwrap_or(1.0);
        let beta = dgrdl_params
            .get("beta")
            .and_then(|v| v.as_real())
            .unwrap_or(1.0);
        let max_iter = dgrdl_params
            .get("max_iter")
            .and_then(|v| v.as_integer())
            .unwrap_or(20) as usize;
        let k_neighbours = dgrdl_params
            .get("k_neighbours")
            .and_then(|v| v.as_integer())
            .unwrap_or(5) as usize;
        let admm_iter = dgrdl_params
            .get("admm_iter")
            .and_then(|v| v.as_integer())
            .unwrap_or(5) as usize;
        let rho = dgrdl_params
            .get("rho")
            .and_then(|v| v.as_real())
            .unwrap_or(1.0);

        DgrdlParams {
            sparsity,
            dict_size,
            alpha,
            beta,
            max_iter,
            k_neighbours,
            admm_iter,
            rho,
        }
    }

    /// Create default DGRDL parameters
    ///
    /// ### Returns
    ///
    /// `DgrdlParams` with default parameter configuration
    fn default() -> Self {
        Self {
            sparsity: 5,
            dict_size: 50,
            alpha: 1.0,
            beta: 1.0,
            max_iter: 20,
            k_neighbours: 5,
            admm_iter: 5,
            rho: 1.0,
        }
    }
}

/// DGRDL algorithm results
///
/// ### Fields
///
/// * `dictionary` - Learned dictionary of size n x k
/// * `coefficients` - Sparse coefficients of size k x m
/// * `context_laplacian` - Sample context Laplacian n × n
/// * `gene_laplacian` - Feature context laplacian of size m x m
#[derive(Debug)]
pub struct DgrdlResults {
    pub dictionary: Mat<f64>,
    pub coefficients: Mat<f64>,
    pub sample_laplacian: Mat<f64>,
    pub feature_laplacian: Mat<f64>,
}

/////////////////
// DGRDL Cache //
/////////////////

/// DGRDL cache for hyperparameter tuning
///
/// ### Fields
///
/// * `data_hash` - Hash of the underlying data
/// * `laplacian_cache` - HashMap of the feature and sample Laplacian matrix
/// * `dictionary_cache` - HashMap of the dictionary cashe with (dict_size, seed)
///   as keys
/// * `distance_matrix` - The distance matrix of the data as a flat vector
#[derive(Debug)]
pub struct DgrdlCache {
    pub data_hash: u64,
    pub laplacian_cache: FxHashMap<usize, (Mat<f64>, Mat<f64>)>,
    pub dictionary_cache: FxHashMap<(usize, usize), Mat<f64>>,
    pub distance_matrix: Option<Vec<f64>>,
}

impl DgrdlCache {
    /// Generate a new instance of the Cache
    ///
    /// ### Returns
    ///
    /// Returns self with empty fields
    pub fn new() -> Self {
        Self {
            data_hash: 0,
            laplacian_cache: FxHashMap::default(),
            dictionary_cache: FxHashMap::default(),
            distance_matrix: None,
        }
    }

    /// Clear cache if data has changed
    ///
    /// Will clear the cash if the Hash changed
    ///
    /// ### Params
    ///
    /// * `data` The data for which to calculate the Hash
    fn validate_data(&mut self, data: &MatRef<f64>) {
        let new_hash = self.compute_data_hash(data);
        if new_hash != self.data_hash {
            self.clear_cache();
            self.data_hash = new_hash;
        }
    }

    /// Compute the hash of the data. This will test just some elements
    ///
    /// ### Params
    ///
    /// * `data` - The data for the DGRDL algorithm
    ///
    /// ### Returns
    ///
    /// Hash usize
    fn compute_data_hash(&self, data: &MatRef<f64>) -> u64 {
        let mut hash = data.nrows() as u64;
        hash = hash.wrapping_mul(31).wrapping_add(data.ncols() as u64);
        if data.nrows() > 0 && data.ncols() > 0 {
            hash = hash.wrapping_mul(31).wrapping_add(data[(0, 0)].to_bits());
            let last_row = data.nrows() - 1;
            let last_col = data.ncols() - 1;
            hash = hash
                .wrapping_mul(31)
                .wrapping_add(data[(last_row, last_col)].to_bits());
        }
        hash
    }

    /// Clears the internal Cache
    fn clear_cache(&mut self) {
        self.laplacian_cache.clear();
        self.dictionary_cache.clear();
        self.distance_matrix = None;
    }
}

///////////////////
// DGRDL objects //
///////////////////

/// Calculate individual components of DGRDL objective function
///
/// Useful for hyperparameter tuning and model selection
///
/// ### Fields
///
/// * `approximation_error` - Calculates the approximation error in form of
///   squared Frobenius norm
/// * `feature_laplacian_objective` - Measures the similiarity of coefficients
///   for the feature Laplacian
/// * `sample_laplacian_objective` - Measures the smoothness of the dictionary
///   in respect to sample Laplacian
/// * `seed` - Random seed for reproducibility
/// * `k_neighbours` - Number of neighbours
/// * `dict_size` - The size of the dictionary
#[derive(Debug, Clone)]
pub struct DgrdlObjectives {
    pub approximation_error: f64,
    pub feature_laplacian_objective: f64,
    pub sample_laplacian_objective: f64,
    pub seed: usize,
    pub k_neighbours: usize,
    pub dict_size: usize,
}

impl DgrdlObjectives {
    /// Calculate the different objectives and reconstruction error for a given DGRDL run
    ///
    /// ### Params
    ///
    /// * `dat` - Input matrix
    /// * `res` - The `DgrdlResults` of that run
    /// * `alpha` - Sample context regularisation weight
    /// * `beta` - Feature context regularisation weight
    ///
    /// ### Returns
    ///
    /// The `DgrdlObjectives` object
    pub fn calculate(
        dat: &MatRef<f64>,
        res: &DgrdlResults,
        alpha: f64,
        beta: f64,
        seed: usize,
        k_neighbours: usize,
        dict_size: usize,
    ) -> Self {
        let approximation_error =
            Self::calculate_approximation_error(dat, &res.dictionary, &res.coefficients);

        let feature_laplacian_objective =
            beta * Self::calculate_trace_xlx(&res.coefficients, &res.feature_laplacian);

        let sample_laplacian_objective =
            alpha * Self::calculate_trace_ddl(&res.dictionary, &res.sample_laplacian);

        Self {
            approximation_error,
            feature_laplacian_objective,
            sample_laplacian_objective,
            seed,
            k_neighbours,
            dict_size,
        }
    }

    /// Calculate the reconstruction error
    ///
    /// ### Params
    ///
    /// * `x` - The input matrix
    /// * `dictionary` - The yielded dictionary matrix from the algorithm
    /// * `coefficients` - The yielded coefficient matrix from the algorithm
    ///
    /// ### Returns
    ///
    /// The Frobenius norm squared. The lower, the better.
    fn calculate_approximation_error(
        x: &MatRef<f64>,
        dictionary: &Mat<f64>,
        coefficients: &Mat<f64>,
    ) -> f64 {
        let reconstruction = dictionary * coefficients;
        let residual = x.as_ref() - &reconstruction;
        residual.norm_l2().powi(2)
    }

    /// Calculate gene Laplacian term: tr(X·L_f·X^T)
    ///
    /// ### Params
    ///
    /// * `coefficients` - The yielded coefficient matrix from the algorithm
    /// * `feature_laplacian` - The Laplacian matrix of the features
    ///
    /// ### Returns
    ///
    /// The trace of that part of the equation. The lower, the better.
    fn calculate_trace_xlx(coefficients: &Mat<f64>, feature_laplacian: &Mat<f64>) -> f64 {
        let xl = coefficients * feature_laplacian;

        let mut trace = 0.0;
        for i in 0..coefficients.nrows() {
            for j in 0..coefficients.ncols() {
                trace += xl[(i, j)] * coefficients[(i, j)];
            }
        }
        trace
    }

    /// Calculate cell Laplacian term: tr(D^T·D·L_s)
    ///
    /// ### Params
    ///
    /// * `dictionary` - The yielded dictionary matrix from the algorithm
    /// * `sample_laplacian` - The Laplacian matrix of the samples
    ///
    /// ### Returns
    ///
    /// The trace of that part of equation. The lower, the better.
    fn calculate_trace_ddl(dictionary: &Mat<f64>, sample_laplacian: &Mat<f64>) -> f64 {
        // D^T·L_s (dict_size × n_samples)
        let dt_l = dictionary.transpose() * sample_laplacian;

        // tr(D^T·L_s·D) = sum_{i,j} (D^T·L_s)[i,j] * D[j,i]
        let mut trace = 0.0;
        for i in 0..dictionary.ncols() {
            // dict_size
            for j in 0..dictionary.nrows() {
                // n_samples
                trace += dt_l[(i, j)] * dictionary[(j, i)];
            }
        }
        trace
    }
}

//////////////////////////
// DGRDL implementation //
//////////////////////////

/// Main DGRDL implementation
///
/// Implements the Dual Graph Regularized Dictionary Learning algorithm which
/// learns a sparse dictionary representation while preserving local geometric
/// structure in both sample and feature spaces through graph regularization.
///
/// ### Fields
///
/// * `params` - Configuration parameters for the algorithm
#[derive(Debug)]
pub struct Dgrdl {
    params: DgrdlParams,
    cache: DgrdlCache,
}

impl Dgrdl {
    /// Create a new DGRDL instance
    ///
    /// ### Params
    ///
    /// * `params` - Configuration parameters for the algorithm
    ///
    /// ### Returns
    ///
    /// New DGRDL instance ready for training
    pub fn new(params: DgrdlParams) -> Self {
        Self {
            params,
            cache: DgrdlCache::new(),
        }
    }

    /// Run DGRDL algorithm on input data (with initial hyper parameters)
    ///
    /// Implements the main optimization loop alternating between sparse coding
    /// and dictionary update steps. The algorithm minimizes:
    /// ||Y - DX||²_F + α·tr(D^T·D·L_s) + β·tr(X·L_f·X^T) + sparsity constraint
    /// where Y is data, D is dictionary, X is coefficients, L_s and L_f are
    /// sample and feature Laplacians respectively.
    ///
    /// ### Params
    ///
    /// * `data` - Input data matrix of size n_samples × n_features
    /// * `seed` - Seed for dictionary initialisation.
    /// * `verbose` - Whether to print iteration progress
    ///
    /// ### Returns
    ///
    /// `DgrdlResults` containing the learned dictionary, coefficients, and Laplacians
    pub fn fit(&mut self, data: &MatRef<f64>, seed: usize, verbose: bool) -> DgrdlResults {
        self.fit_with_params(data, self.params.clone(), seed, verbose)
    }

    /// Run DGRDL algorithm on input data (with provided hyperparameters)
    ///
    /// Implements the main optimization loop alternating between sparse coding
    /// and dictionary update steps. The algorithm minimizes:
    /// ||Y - DX||²_F + α·tr(D^T·D·L_s) + β·tr(X·L_f·X^T) + sparsity constraint
    /// where Y is data, D is dictionary, X is coefficients, L_s and L_f are
    /// sample and feature Laplacians respectively.
    ///
    /// ### Params
    ///
    /// * `data` - Input data matrix of size n_samples × n_features
    /// * `params` - The `DgrdlParams` containing the parameters you wish to use.
    /// * `seed` - Seed for dictionary initialisation.
    /// * `verbose` - Whether to print iteration progress
    ///
    /// ### Returns
    ///
    /// `DgrdlResults` containing the learned dictionary, coefficients, and Laplacians
    pub fn fit_with_params(
        &mut self,
        data: &MatRef<f64>,
        params: DgrdlParams,
        seed: usize,
        verbose: bool,
    ) -> DgrdlResults {
        let (n_samples, n_features) = data.shape();

        if verbose {
            println!(
                "Starting DGRDL run with {:?} samples and {:?} features.",
                n_samples, n_features
            )
        }

        let start_total = Instant::now();

        let start_dictionary_gen = Instant::now();
        let mut dictionary = self.get_dictionary(data, params.dict_size, seed);
        let end_dictionary_gen = start_dictionary_gen.elapsed();

        if verbose {
            println!(
                "DGRDL dictionary obtained in {:.2?} (cached: {}).",
                end_dictionary_gen,
                self.cache
                    .dictionary_cache
                    .contains_key(&(self.params.dict_size, seed))
            );
        }

        let start_laplacian_gen = Instant::now();
        let (feature_laplacian, sample_laplacian) = self.get_laplacians(data, params.k_neighbours);
        let end_laplacian_gen = start_laplacian_gen.elapsed();

        if verbose {
            println!(
                "DGRDL graph laplacians obtained in {:.2?} (cached: {}).",
                end_laplacian_gen,
                self.cache
                    .laplacian_cache
                    .contains_key(&params.k_neighbours)
            );
        }

        let mut coefficients: Mat<f64> = Mat::zeros(params.dict_size, n_features);

        for iter in 0..self.params.max_iter {
            let start_iter = Instant::now();

            // Sparse coding step
            coefficients = grsc_admm(
                &dictionary,
                data,
                &feature_laplacian,
                params.sparsity,
                params.beta,
                params.admm_iter,
                params.rho,
            );

            // Dictionary update step
            dictionary = update_dictionary(data, &coefficients, &sample_laplacian, params.alpha);

            let end_iter = start_iter.elapsed();

            if verbose {
                println!(
                    " DGRDL iteration {}/{} in {:.2?}.",
                    iter + 1,
                    params.max_iter,
                    end_iter
                );
            }
        }

        let end_time = start_total.elapsed();

        if verbose {
            println!(
                "Total time elapsed for fitting the DGRDL run: {:.2?}.",
                end_time
            )
        }

        DgrdlResults {
            dictionary,
            coefficients,
            sample_laplacian,
            feature_laplacian,
        }
    }

    pub fn grid_search(
        &mut self,
        data: &MatRef<f64>,
        dict_sizes: &[usize],
        k_neighbours_iters: &[usize],
        seeds: &[usize],
        verbose: bool,
    ) -> Vec<DgrdlObjectives> {
        let (n_samples, n_features) = data.shape();

        let n_dict = dict_sizes.len();
        let n_k_neighbours = k_neighbours_iters.len();
        let n_seeds = k_neighbours_iters.len();

        let total_permutations = n_dict * n_k_neighbours * n_seeds;

        if verbose {
            println!(
                "Starting DGRDL hyperparameter grid search with {:?} samples and {:?} features and {:?} permutations.",
                n_samples, n_features, total_permutations
            )
        }

        let start_total = Instant::now();
        let mut iteration: usize = 1;

        let mut grid_search_res: Vec<DgrdlObjectives> = Vec::with_capacity(total_permutations);

        for seed in seeds {
            for dict_size in dict_sizes {
                for &k_neighbours in k_neighbours_iters {
                    let iter_total = Instant::now();

                    let mut params = self.params.clone();
                    params.dict_size = *dict_size;
                    params.k_neighbours = k_neighbours;

                    if verbose {
                        println!(
                            " Iter ({}|{}) - seed = {} | dict_size = {} | k_neighbours = {}",
                            iteration, total_permutations, seed, dict_size, k_neighbours
                        );
                    }

                    let res = self.fit_with_params(data, params.clone(), *seed, false);

                    let metrics = DgrdlObjectives::calculate(
                        data,
                        &res,
                        params.alpha,
                        params.beta,
                        *seed,
                        k_neighbours,
                        *dict_size,
                    );

                    grid_search_res.push(metrics);

                    let iter_time = iter_total.elapsed();

                    if verbose {
                        println!(
                            " Iter ({}|{}) - finalised in {:.2?}",
                            iteration, total_permutations, iter_time
                        );
                    }

                    iteration += 1;
                }
            }
        }

        let total_time = start_total.elapsed();

        if verbose {
            println!("Finished hyperparameter grid search in {:.2?}", total_time);
        }

        grid_search_res
    }

    /// Get or compute dictionary for given dict_size and seed
    ///
    /// ### Params
    ///
    /// * `data` - Input data into the algorithm
    /// * `dict_size` - The size of the dictionary
    /// * `seed` - Random seed for the initialisation
    ///
    /// ### Returns
    ///
    /// The initial dictionary matrix
    fn get_dictionary(&mut self, data: &MatRef<f64>, dict_size: usize, seed: usize) -> Mat<f64> {
        self.cache.validate_data(data);

        let key = (dict_size, seed);
        if let Some(dictionary) = self.cache.dictionary_cache.get(&key) {
            return dictionary.clone();
        }

        // Ensure distance matrix is computed (reused across dictionary initializations)
        if self.cache.distance_matrix.is_none() {
            self.cache.distance_matrix = Some(self.precompute_distances(data));
        }

        // Compute new dictionary using cached distances
        let dictionary = self.initialise_dictionary(data, dict_size, seed);

        // Cache the result
        self.cache.dictionary_cache.insert(key, dictionary.clone());

        dictionary
    }

    /// Initialise dictionary using k-medoids clustering
    ///
    /// Selects initial dictionary atoms by finding representative samples
    /// (medoids) that minimise the total distance to all other samples.
    /// This provides a good initialisation that captures data diversity. Will
    /// leverage the pre-computed distances.
    ///
    /// ### Params
    ///
    /// * `data` - Input data matrix
    /// * `dict_size` - The size of the dictionary.
    /// * `seed` - Seed for the random initilisation of the first medoid.
    ///
    /// ### Returns
    ///
    /// Initialised dictionary matrix with normalized columns
    fn initialise_dictionary(&self, data: &MatRef<f64>, dict_size: usize, seed: usize) -> Mat<f64> {
        let mut rng = StdRng::seed_from_u64(seed as u64);

        let k = dict_size;
        let (n_samples, n_features) = data.shape();

        let mut selected = Vec::with_capacity(k);

        // Choose first center randomly
        selected.push(rng.random_range(0..n_features));

        // Precompute all pairwise distances (only upper triangle)
        let distances = self.cache.distance_matrix.as_ref().unwrap();

        for _ in 1..k {
            // Find distances to nearest selected center
            let min_distances: Vec<f64> = (0..n_features)
                .into_par_iter()
                .map(|i| {
                    if selected.contains(&i) {
                        0.0 // Already selected
                    } else {
                        selected
                            .iter()
                            .map(|&j| get_distance(distances, i, j, n_features))
                            .fold(f64::INFINITY, f64::min)
                    }
                })
                .collect();

            // Choose next center with probability proportional to squared distance
            let total_weight: f64 = min_distances.iter().map(|d| d * d).sum();

            if total_weight > 0.0 {
                let mut target = rng.random::<f64>() * total_weight;
                let mut chosen = 0;

                for (i, &dist) in min_distances.iter().enumerate() {
                    target -= dist * dist;
                    if target <= 0.0 || i == n_features - 1 {
                        chosen = i;
                        break;
                    }
                }

                if !selected.contains(&chosen) {
                    selected.push(chosen);
                }
            } else {
                // Fallback: choose randomly from unselected
                let unselected: Vec<usize> =
                    (0..n_features).filter(|i| !selected.contains(i)).collect();
                if !unselected.is_empty() {
                    selected.push(unselected[rng.random_range(0..unselected.len())]);
                }
            }
        }

        // Build dictionary
        let mut dictionary = Mat::zeros(n_samples, k);
        for (dict_idx, &feature_idx) in selected.iter().enumerate() {
            let col = data.col(feature_idx);
            let norm = col.norm_l2();
            if norm > 1e-12 {
                for i in 0..n_samples {
                    dictionary[(i, dict_idx)] = col[i] / norm;
                }
            }
        }
        dictionary
    }

    /// Helper function to pre-calculate all the distances
    ///
    /// This function will calculate the cosine distance between all the
    /// columns in a fast, sparse way
    ///
    /// ### Params
    ///
    /// * `data` - The input data
    ///
    /// ### Returns
    ///
    /// The flattened distances
    fn precompute_distances(&self, data: &MatRef<f64>) -> Vec<f64> {
        let n_features = data.ncols();
        let n_pairs = n_features * (n_features - 1) / 2;

        (0..n_pairs)
            .into_par_iter()
            .map(|pair_idx| {
                let (i, j) = triangle_to_indices(pair_idx, n_features);
                column_distance(data.col(i), data.col(j))
            })
            .collect()
    }

    /// Pull out the Laplacians out of the cache or recompute for given k_neighbours
    ///
    /// ### Params
    ///
    /// * `data` - Input data
    /// * `k_neighbours` - Number of neighbours in the KNN graph
    ///
    /// ### Returns
    ///
    /// Tuple of `(feature_laplacian, sample_laplacian)`
    fn get_laplacians(&mut self, data: &MatRef<f64>, k_neighbours: usize) -> (Mat<f64>, Mat<f64>) {
        self.cache.validate_data(data);

        if let Some((feature_lap, sample_lap)) = self.cache.laplacian_cache.get(&k_neighbours) {
            return (feature_lap.clone(), sample_lap.clone());
        }

        // Compute new Laplacians
        let feature_laplacian = get_dgrdl_laplacian(&data.as_ref(), k_neighbours, true);
        let sample_laplacian = get_dgrdl_laplacian(&data.as_ref(), k_neighbours, false);

        // Cache the results
        self.cache.laplacian_cache.insert(
            k_neighbours,
            (feature_laplacian.clone(), sample_laplacian.clone()),
        );

        (feature_laplacian, sample_laplacian)
    }
}

/////////////
// Helpers //
/////////////

/// Get the upper triangle indices as a pair for rapid distance calculations
///
/// Stores the indices of a 3 x 3 matrix for example in a pattern of
/// `(0,1)`, `(0,2)`, `(0,3)`, `(1,2)`, `(1,3)`, `(2,3)` (assuming a linear
/// matrix).
///
/// ### Params
///
/// * `pair_idx` - Linear index in compressed storage, i.e., `0 to n (n - 1) / 2 - 1)`
/// * `n` - Shape of the original column
///
/// ### Returns
///
/// Tuple of `(i, j)`
fn triangle_to_indices(linear_idx: usize, n: usize) -> (usize, usize) {
    let mut idx = linear_idx;
    let mut i = 0;

    while idx >= n - i - 1 {
        idx -= n - i - 1;
        i += 1;
    }

    (i, i + 1 + idx)
}

/// Retrieve distance from compressed upper-triangle storage
///
/// ### Params
///
/// * `distances` - Pre-computed distance array in compressed format
/// * `i` - First matrix index
/// * `j` - Second matrix index  
/// * `n` - Original matrix size
///
/// ### Returns
///
/// Distance value between the two values
fn get_distance(distances: &[f64], i: usize, j: usize, n: usize) -> f64 {
    if i == j {
        return 0.0;
    }

    let (min_idx, max_idx) = if i < j { (i, j) } else { (j, i) };
    let linear_idx = (min_idx * (2 * n - min_idx - 1)) / 2 + (max_idx - min_idx - 1);
    distances[linear_idx]
}

/// Create the Laplacian matrix for the features or samples
///
/// ### Fields
///
/// * `data` - The original data matrix with rows = samples and columns =
///   features
/// * `k` - Number of neighbours for the KNN graph
/// * `features` - If `true` generated the Laplacian for the features, otherwise
///   for the samples
///
/// ### Returns
///
/// The Laplacian matrix L = D - A where D is the degree matrix and A is the
/// adjacency matrix
fn get_dgrdl_laplacian(data: &MatRef<f64>, k: usize, features: bool) -> Mat<f64> {
    let cosine_sim = if features {
        column_cosine(data)
    } else {
        column_cosine(&data.transpose())
    };

    let knn_adjacency = get_knn_graph_adj(&cosine_sim.as_ref(), k);

    adjacency_to_laplacian(&knn_adjacency.as_ref())
}

/// Function to do the ADMM on the GRSC
///
/// Solves the optimization problem:
///
/// min_X ||Y - DX||²_F + β·tr(X·L_f·X^T) + λ||X||_0
///
/// subject to sparsity constraint using alternating direction method of
/// multipliers (ADMM).
///
/// ### Params
///
/// * `dictionary` - Current dictionary matrix D
/// * `data` - Input data matrix Y
/// * `feature_laplacian` - Feature graph Laplacian L_f
/// * `sparsity` - Maximum number of non-zero coefficients per column
/// * `beta` - Graph regularization weight
/// * `max_iter` - Maximum ADMM iterations
/// * `rho` - ADMM penalty parameter
///
/// ### Returns
///
/// Sparse coefficient matrix X satisfying the sparsity constraint
pub fn grsc_admm(
    dictionary: &Mat<f64>,
    data: &MatRef<f64>,
    feature_laplacian: &Mat<f64>,
    sparsity: usize,
    beta: f64,
    max_iter: usize,
    rho: f64,
) -> Mat<f64> {
    let k = dictionary.ncols();
    let m = data.ncols();

    let gram = dictionary.transpose() * dictionary;
    let dt_y = dictionary.transpose() * data;

    // Initialize ADMM variables
    #[allow(unused_assignments)]
    let mut x: Mat<f64> = Mat::zeros(k, m);
    let mut z: Mat<f64> = Mat::zeros(k, m);
    let mut u: Mat<f64> = Mat::zeros(k, m);

    for _ in 0..max_iter {
        // X update: solve (G + ρI)X + βXL = D^T Y + ρ(Z - U)
        x = admm_solve_x(&gram, &dt_y, &z, &u, feature_laplacian, beta, rho);

        // Z update: sparse projection
        z = sparse_projection(&(&x + &u), sparsity);

        // U update: dual variable
        u = &u + &x - &z;
    }

    z
}

/// Solve the X-update step in ADMM
///
/// Solves the Sylvester equation: AX + XB = C
/// where A = G + ρI, B = βL_f, C = D^T Y + ρ(Z - U)
///
/// ### Params
///
/// * `gram` - Gram matrix G = D^T D
/// * `dt_y` - Data projection D^T Y
/// * `z` - Current Z variable from ADMM
/// * `u` - Current dual variable from ADMM
/// * `feature_laplacian` - Feature Laplacian L_f
/// * `beta` - Graph regularization weight
/// * `rho` - ADMM penalty parameter
///
/// ### Returns
///
/// Updated X matrix solving the regularized least squares problem
fn admm_solve_x(
    gram: &Mat<f64>,
    dt_y: &Mat<f64>,
    z: &Mat<f64>,
    u: &Mat<f64>,
    feature_laplacian: &Mat<f64>,
    beta: f64,
    rho: f64,
) -> Mat<f64> {
    let k = gram.nrows();

    let identity_mat: Mat<f64> = Mat::identity(k, k);
    let a: Mat<f64> = gram + &(rho * identity_mat);
    let b = beta * feature_laplacian;
    let c = dt_y + &(rho * (z - u));

    sylvester_solver(&a.as_ref(), &b.as_ref(), &c.as_ref())
}

/// Apply sparsity constraint via hard thresholding
///
/// For each column, keeps only the largest (in absolute value) 'sparsity'
/// number of coefficients and sets the rest to zero.
///
/// ### Params
///
/// * `x` - Input coefficient matrix
/// * `sparsity` - Number of non-zero coefficients to keep per column
///
/// ### Returns
///
/// Sparse matrix with at most 'sparsity' non-zero entries per column
fn sparse_projection(x: &Mat<f64>, sparsity: usize) -> Mat<f64> {
    let (k, m) = x.shape();
    let mut result = Mat::zeros(k, m);

    // Process columns in parallel
    let columns: Vec<Vec<f64>> = (0..m)
        .into_par_iter()
        .map(|j| {
            let x_col = x.col(j);
            let mut indexed_vals: Vec<(usize, f64)> =
                x_col.iter().enumerate().map(|(i, &val)| (i, val)).collect();

            // Partial sort to find top sparsity elements
            let take_count = sparsity.min(k);
            if take_count < indexed_vals.len() {
                indexed_vals.select_nth_unstable_by(take_count, |a, b| {
                    b.1.abs().partial_cmp(&a.1.abs()).unwrap()
                });
            }

            // Create column result
            let mut col_result = vec![0.0; k];
            for &(idx, val) in indexed_vals.iter().take(sparsity) {
                col_result[idx] = val;
            }
            col_result
        })
        .collect();

    // Copy results back to matrix
    for (j, col_data) in columns.iter().enumerate() {
        for (i, &val) in col_data.iter().enumerate() {
            result[(i, j)] = val;
        }
    }

    result
}

/// Compute Euclidean distance between two columns
///
/// ### Params
///
/// * `col_i` - First column vector
/// * `col_j` - Second column vector
///
/// ### Returns
///
/// Euclidean distance between the two columns
fn column_distance(col_i: ColRef<f64>, col_j: ColRef<f64>) -> f64 {
    col_i
        .iter()
        .zip(col_j.iter())
        .map(|(a, b)| (a - b).powi(2))
        .sum::<f64>()
        .sqrt()
}

/// Update the entire dictionary
///
/// Updates all dictionary atoms sequentially using graph regularisation.
/// Each atom is updated to minimise the reconstruction error while
/// maintaining smoothness according to the sample graph structure.
///
/// ### Params
///
/// * `data` - Input data matrix Y
/// * `coefficients` - Current coefficient matrix X
/// * `sample_laplacian` - Sample graph Laplacian L_s
/// * `alpha` - Sample regularization weight
/// * `beta` - Feature regularization weight
///
/// ### Returns
///
/// Updated dictionary matrix with normalised columns
fn update_dictionary(
    data: &MatRef<f64>,
    coefficients: &Mat<f64>,
    sample_laplacian: &Mat<f64>,
    alpha: f64,
) -> Mat<f64> {
    let (n_contexts, _) = data.shape();
    let k = coefficients.nrows();

    // Precompute common matrices
    let identity: Mat<f64> = Mat::identity(n_contexts, n_contexts);
    let reg_term = alpha * sample_laplacian;

    // Process atoms in parallel, collect results
    let atom_columns: Vec<Vec<f64>> = (0..k)
        .into_par_iter()
        .map(|atom_idx| {
            let x_j = coefficients.row(atom_idx);

            // Find active signals
            let active_signals: Vec<(usize, f64)> = x_j
                .iter()
                .enumerate()
                .filter_map(|(signal_idx, &coeff)| {
                    if coeff.abs() > 1e-12 {
                        Some((signal_idx, coeff))
                    } else {
                        None
                    }
                })
                .collect();

            if active_signals.is_empty() {
                return vec![0.0; n_contexts];
            }

            // Compute right-hand side
            let mut rhs = Mat::zeros(n_contexts, 1);
            for &(signal_idx, coeff) in &active_signals {
                let signal = data.col(signal_idx);
                for i in 0..n_contexts {
                    rhs[(i, 0)] += signal[i] * coeff;
                }
            }

            // Solve system
            let x_j_norm_sq: f64 = active_signals.iter().map(|(_, coeff)| coeff * coeff).sum();

            if x_j_norm_sq > 1e-12 {
                let system_matrix = &(x_j_norm_sq * &identity) + &reg_term;
                let lu = PartialPivLu::new(system_matrix.as_ref());
                let solution = lu.solve(&rhs);

                // Normalize
                let norm = solution.norm_l2();
                if norm > 1e-12 {
                    (0..n_contexts).map(|i| solution[(i, 0)] / norm).collect()
                } else {
                    vec![0.0; n_contexts]
                }
            } else {
                vec![0.0; n_contexts]
            }
        })
        .collect();

    // Build dictionary from results
    let mut dictionary = Mat::zeros(n_contexts, k);
    for (atom_idx, atom_col) in atom_columns.iter().enumerate() {
        for (i, &val) in atom_col.iter().enumerate() {
            dictionary[(i, atom_idx)] = val;
        }
    }

    dictionary
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::prelude::*;

    #[test]
    fn test_known_dictionary_recovery() {
        let mut rng = StdRng::seed_from_u64(42);

        // Create ground truth dictionary
        let n_samples = 15;
        let n_features = 30;
        let true_dict_size = 4;
        let sparsity = 2;

        let mut true_dict = Mat::zeros(n_samples, true_dict_size);

        // Create 4 simple distinct patterns
        // Pattern 1: Constant high
        for i in 0..n_samples {
            true_dict[(i, 0)] = 1.0;
        }

        // Pattern 2: Linear
        for i in 0..n_samples {
            true_dict[(i, 1)] = (i as f64) / (n_samples as f64);
        }

        // Pattern 3: Step function
        for i in 0..n_samples {
            true_dict[(i, 2)] = if i < n_samples / 2 { 1.0 } else { -1.0 };
        }

        // Pattern 4: Alternating
        for i in 0..n_samples {
            true_dict[(i, 3)] = if i % 2 == 0 { 1.0 } else { -1.0 };
        }

        // Normalise columns
        for j in 0..true_dict_size {
            let col = true_dict.col(j);
            let norm = col.norm_l2();
            for i in 0..n_samples {
                true_dict[(i, j)] /= norm;
            }
        }

        // Generate sparse coefficients
        let mut true_coeffs = Mat::zeros(true_dict_size, n_features);
        for j in 0..n_features {
            // Randomly select which atoms to use
            let mut indices: Vec<usize> = (0..true_dict_size).collect();
            indices.shuffle(&mut rng);

            for &idx in indices.iter().take(sparsity) {
                true_coeffs[(idx, j)] = rng.random_range(0.5..2.0); // Positive coefficients only
            }
        }

        // Generate synthetic data: Y = D * X + small noise
        let mut synthetic_data = &true_dict * &true_coeffs;

        // Add minimal noise
        for i in 0..n_samples {
            for j in 0..n_features {
                synthetic_data[(i, j)] += rng.random_range(-0.05..0.05);
            }
        }

        // Run DGRDL with fast settings
        // Weak regularisation, less iterations
        let params = DgrdlParams {
            sparsity,
            dict_size: true_dict_size,
            alpha: 0.25,
            beta: 0.25,
            max_iter: 8,
            k_neighbours: 3,
            admm_iter: 5,
            rho: 1.0,
        };

        let mut dgrdl = Dgrdl::new(params);
        let result = dgrdl.fit(&synthetic_data.as_ref(), 42, true);

        // Test reconstruction quality
        let reconstruction = &result.dictionary * &result.coefficients;
        let reconstruction_error =
            (&synthetic_data - &reconstruction).norm_l2() / synthetic_data.norm_l2();

        println!("Reconstruction error: {:.6}", reconstruction_error);
        assert!(
            reconstruction_error < 0.3,
            "Poor reconstruction: {}",
            reconstruction_error
        ); // Relaxed threshold

        // Test sparsity
        let avg_sparsity: f64 = (0..n_features)
            .map(|j| {
                result
                    .coefficients
                    .col(j)
                    .iter()
                    .filter(|&&x| x.abs() > 1e-6)
                    .count() as f64
            })
            .sum::<f64>()
            / n_features as f64;

        println!(
            "Average sparsity: {:.2} (target: {})",
            avg_sparsity, sparsity
        );
        assert!(
            avg_sparsity <= sparsity as f64 + 1.0,
            "Sparsity constraint violated"
        );
    }

    /// Test 2: Block-structured data (simulates gene modules)
    #[test]
    fn test_block_structure() {
        let n_samples = 20;
        let n_features = 40;
        let n_blocks = 4;
        let features_per_block = n_features / n_blocks;

        let mut data = Mat::zeros(n_samples, n_features);

        // Create 4 blocks with different patterns
        for block in 0..n_blocks {
            let start_sample = block * (n_samples / n_blocks);
            let end_sample = (block + 1) * (n_samples / n_blocks);
            let start_feature = block * features_per_block;
            let end_feature = (block + 1) * features_per_block;

            // Fill block with correlated values
            let base_value = (block as f64 + 1.0) * 0.5;
            for i in start_sample..end_sample {
                for j in start_feature..end_feature {
                    data[(i, j)] = base_value + 0.1 * ((i + j) as f64).sin();
                }
            }
        }

        // Add minimal noise
        let mut rng = StdRng::seed_from_u64(123);
        for i in 0..n_samples {
            for j in 0..n_features {
                data[(i, j)] += rng.random_range(-0.05..0.05);
            }
        }

        let params = DgrdlParams {
            sparsity: 2,
            dict_size: 6,
            alpha: 0.5,
            beta: 0.5,
            max_iter: 8,
            k_neighbours: 3,
            admm_iter: 5,
            rho: 1.0,
        };

        let mut dgrdl = Dgrdl::new(params);
        let result = dgrdl.fit(&data.as_ref(), 42, true);

        // Check that we can reconstruct the block structure
        let reconstruction = &result.dictionary * &result.coefficients;
        let block_error = (&data - &reconstruction).norm_l2() / data.norm_l2();

        println!("Block reconstruction error: {:.6}", block_error);
        assert!(block_error < 0.5, "Failed to capture block structure");
    }

    /// Test 3: Test something that represents biology
    #[test]
    fn test_gene_modules() {
        let n_cell_lines = 24;
        let n_genes = 60;

        let mut data = Mat::zeros(n_cell_lines, n_genes);

        // Create modules representing biological functions
        let modules = [
            ("DNA_repair", 0..15),  // Module 1
            ("Cell_cycle", 15..30), // Module 2
            ("Metabolism", 30..45), // Module 3
            ("Stress", 45..60),     // Module 4
        ];

        let mut rng = StdRng::seed_from_u64(456);

        // Simulate module activities across cell lines
        for (module_idx, (name, gene_range)) in modules.iter().enumerate() {
            println!("Creating module: {}", name);

            // Each module is active in different cell lines
            let active_cells = match module_idx {
                0 => 0..8,   // Have module 1 active
                1 => 6..14,  // Have module 2 active
                2 => 12..20, // Have module 3 active
                3 => 16..24, // Have module 4 active
                _ => 0..0,
            };

            // Set module activity
            for cell in active_cells {
                let base_activity = 1.0 + rng.random_range(-0.2..0.2);
                for gene in gene_range.clone() {
                    // Genes in same module have correlated activities
                    data[(cell, gene)] = base_activity + rng.random_range(-0.1..0.1);
                }
            }
        }

        // Add background noise
        for i in 0..n_cell_lines {
            for j in 0..n_genes {
                data[(i, j)] += rng.random_range(-0.05..0.05);
            }
        }

        let params = DgrdlParams {
            sparsity: 3,
            dict_size: 8,
            alpha: 0.5,
            beta: 0.5,
            max_iter: 10,
            k_neighbours: 3,
            admm_iter: 5,
            rho: 1.0,
        };

        let mut dgrdl = Dgrdl::new(params);
        let result = dgrdl.fit(&data.as_ref(), 42, true);

        // Validate results
        let reconstruction = &result.dictionary * &result.coefficients;
        let module_error = (&data - &reconstruction).norm_l2() / data.norm_l2();

        println!("Gene module reconstruction error: {:.6}", module_error);

        // Quick check of module similarity (simplified)
        let mut total_within_sim = 0.0;
        let mut count = 0;

        for (_, gene_range) in modules.iter() {
            for gene1 in gene_range.clone().take(3) {
                // Only check first 3 genes per module
                for gene2 in gene_range.clone().skip(1).take(2) {
                    // Check against 2 others
                    if gene1 < gene2 {
                        let coeff1 = result.coefficients.col(gene1);
                        let coeff2 = result.coefficients.col(gene2);

                        let dot_product: f64 =
                            coeff1.iter().zip(coeff2.iter()).map(|(a, b)| a * b).sum();
                        let norm1 = coeff1.norm_l2();
                        let norm2 = coeff2.norm_l2();

                        if norm1 > 1e-6 && norm2 > 1e-6 {
                            total_within_sim += dot_product / (norm1 * norm2);
                            count += 1;
                        }
                    }
                }
            }
        }

        if count > 0 {
            let avg_sim = total_within_sim / count as f64;
            println!("Average within-module similarity: {:.3}", avg_sim);
        }

        assert!(
            module_error < 0.6,
            "Failed to capture gene module structure"
        );
    }

    /// Test 4: Test graph regularisation
    ///
    /// This test is slow as it has to run over three different settings
    #[test]
    fn test_regularisation_effects() {
        // Create data where graph regularization should help
        let n_samples = 30;
        let n_features = 60;

        let mut data = Mat::zeros(n_samples, n_features);
        let mut rng = StdRng::seed_from_u64(789);

        // Create smooth patterns that should benefit from graph regularisation
        for i in 0..n_samples {
            for j in 0..n_features {
                // Smooth function of both sample and feature indices
                let sample_effect =
                    (i as f64 / n_samples as f64 * 2.0 * std::f64::consts::PI).sin();
                let feature_effect = (j as f64 / n_features as f64 * std::f64::consts::PI).cos();
                data[(i, j)] = sample_effect * feature_effect + rng.random_range(-0.1..0.1);
            }
        }

        // Test with different regularisation strengths
        let test_cases = vec![
            ("No regularization", 0.0, 0.0),
            ("Weak regularization", 0.1, 0.1),
            ("Strong regularization", 1.0, 1.0),
        ];

        for (name, alpha, beta) in test_cases {
            println!("\nTesting: {}", name);

            let params = DgrdlParams {
                sparsity: 3,
                dict_size: 8,
                alpha,
                beta,
                max_iter: 10,
                k_neighbours: 5,
                admm_iter: 5,
                rho: 1.0,
            };

            let mut dgrdl = Dgrdl::new(params);
            let result = dgrdl.fit(&data.as_ref(), 42, true);

            // Check reconstruction
            let reconstruction = &result.dictionary * &result.coefficients;
            let error = (&data - &reconstruction).norm_l2() / data.norm_l2();

            // Check sparsity
            let avg_sparsity: f64 = (0..n_features)
                .map(|j| {
                    result
                        .coefficients
                        .col(j)
                        .iter()
                        .filter(|&&x| x.abs() > 1e-6)
                        .count() as f64
                })
                .sum::<f64>()
                / n_features as f64;

            println!("  Reconstruction error: {:.6}", error);
            println!("  Average sparsity: {:.2}", avg_sparsity);

            assert!(error < 0.5, "Poor reconstruction with {}", name);
            assert!(
                avg_sparsity <= 4.0,
                "Sparsity constraint violated with {}",
                name
            );
        }
    }
}
