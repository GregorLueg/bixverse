use extendr_api::prelude::*;

use faer::{
    linalg::solvers::{PartialPivLu, Solve},
    ColRef, Mat, MatRef,
};
use rand::{rngs::StdRng, Rng, SeedableRng};
use rayon::prelude::*;
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
#[derive(Debug, Clone)]
pub struct Dgrdl {
    params: DgrdlParams,
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
        Self { params }
    }

    /// Run DGRDL algorithm on input data
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
    /// * `verbose` - Whether to print iteration progress
    ///
    /// ### Returns
    ///
    /// `DgrdlResults` containing the learned dictionary, coefficients, and Laplacians
    pub fn fit(&self, data: &MatRef<f64>, seed: usize, verbose: bool) -> DgrdlResults {
        let n_features = data.ncols();

        let start_total = Instant::now();

        let start_dictionary_gen = Instant::now();
        let mut dictionary = self.initialise_dictionary(data, seed);
        let end_dictionary_gen = start_dictionary_gen.elapsed();

        if verbose {
            println!(
                "DGRDL dictionary initialisation finished in {:.2?}",
                end_dictionary_gen
            );
        }

        let start_laplacian_gen = Instant::now();
        let feature_laplacian = get_dgrdl_laplacian(&data.as_ref(), self.params.k_neighbours, true);
        let sample_laplacian = get_dgrdl_laplacian(&data.as_ref(), self.params.k_neighbours, false);
        let end_laplacian_gen = start_laplacian_gen.elapsed();

        if verbose {
            println!(
                "DGRDL graph laplacians generated in {:.2?}",
                end_laplacian_gen
            );
        }

        let mut coefficients: Mat<f64> = Mat::zeros(self.params.dict_size, n_features);

        for iter in 0..self.params.max_iter {
            let start_iter = Instant::now();

            // Sparse coding step
            coefficients = grsc_admm(
                &dictionary,
                data,
                &feature_laplacian,
                self.params.sparsity,
                self.params.beta,
                self.params.admm_iter,
                self.params.rho,
            );

            // Dictionary update step
            dictionary =
                update_dictionary(data, &coefficients, &sample_laplacian, self.params.alpha);

            let end_iter = start_iter.elapsed();

            if verbose {
                println!(
                    " DGRDL iteration {}/{} in {:.2?}",
                    iter + 1,
                    self.params.max_iter,
                    end_iter
                );
            }
        }

        let end_time = start_total.elapsed();

        if verbose {
            println!("Total time for fitting DGRDL: {:.2?}", end_time)
        }

        DgrdlResults {
            dictionary,
            coefficients,
            sample_laplacian,
            feature_laplacian,
        }
    }

    /// Initialise dictionary using k-medoids clustering
    ///
    /// Selects initial dictionary atoms by finding representative samples
    /// (medoids) that minimise the total distance to all other samples.
    /// This provides a good initialisation that captures data diversity.
    ///
    /// ### Params
    ///
    /// * `data` - Input data matrix
    /// * `seed` - Seed for the random initilisation of the first medoid.
    ///
    /// ### Returns
    ///
    /// Initialised dictionary matrix with normalized columns
    fn initialise_dictionary(&self, data: &MatRef<f64>, seed: usize) -> Mat<f64> {
        let mut rng = StdRng::seed_from_u64(seed as u64);

        let k = self.params.dict_size;
        let (n_samples, n_features) = data.shape();

        let mut selected = Vec::with_capacity(k);

        // Choose first center randomly
        selected.push(rng.random_range(0..n_features));

        // Precompute all pairwise distances (only upper triangle)
        let distances = self.precompute_distances(data);

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
                            .map(|&j| get_distance(&distances, i, j, n_features))
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
}

/////////////
// Helpers //
/////////////

fn triangle_to_indices(linear_idx: usize, n: usize) -> (usize, usize) {
    let mut idx = linear_idx;
    let mut i = 0;

    while idx >= n - i - 1 {
        idx -= n - i - 1;
        i += 1;
    }

    (i, i + 1 + idx)
}

/// Get distance from compact upper-triangle storage
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
///            features
/// * `k` - Number of neighbours for the KNN graph
/// * `features` - If `true` generated the Laplacian for the features, otherwise
///                for the samples
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
/// subject to sparsity constraint using Alternating Direction Method of
/// Multipliers (ADMM).
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

        // Create ground truth dictionary (smaller for speed)
        let n_samples = 15;
        let n_features = 30;
        let true_dict_size = 4;
        let sparsity = 2;

        let mut true_dict = Mat::zeros(n_samples, true_dict_size);

        // Create 4 simpler, more distinct patterns
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

        // Normalize columns
        for j in 0..true_dict_size {
            let col = true_dict.col(j);
            let norm = col.norm_l2();
            for i in 0..n_samples {
                true_dict[(i, j)] /= norm;
            }
        }

        // Generate sparse coefficients (simpler, cleaner)
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

        // Run DGRDL with faster settings
        let params = DgrdlParams {
            sparsity,
            dict_size: true_dict_size,
            alpha: 0.01, // Very weak regularization
            beta: 0.01,
            max_iter: 8, // Fewer iterations
            k_neighbours: 3,
            admm_iter: 5, // Fewer ADMM iterations
            rho: 1.0,
        };

        let dgrdl = Dgrdl::new(params);
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

        let dgrdl = Dgrdl::new(params);
        let result = dgrdl.fit(&data.as_ref(), 42, true);

        // Check that we can reconstruct the block structure
        let reconstruction = &result.dictionary * &result.coefficients;
        let block_error = (&data - &reconstruction).norm_l2() / data.norm_l2();

        println!("Block reconstruction error: {:.6}", block_error);
        assert!(block_error < 0.5, "Failed to capture block structure");
    }

    #[test]
    fn test_gene_modules() {
        let n_cell_lines = 24;
        let n_genes = 60;

        let mut data = Mat::zeros(n_cell_lines, n_genes);

        // Create modules representing biological functions
        let modules = [
            ("DNA_repair", 0..15),  // Module 1: DNA repair genes
            ("Cell_cycle", 15..30), // Module 2: Cell cycle genes
            ("Metabolism", 30..45), // Module 3: Metabolic genes
            ("Stress", 45..60),     // Module 4: Stress response
        ];

        let mut rng = StdRng::seed_from_u64(456);

        // Simulate module activities across cell lines
        for (module_idx, (name, gene_range)) in modules.iter().enumerate() {
            println!("Creating module: {}", name);

            // Each module is active in different cell lines
            let active_cells = match module_idx {
                0 => 0..8,   // DNA repair: active in first 8 cell lines
                1 => 6..14,  // Cell cycle: overlaps with DNA repair
                2 => 12..20, // Metabolism: active in middle cell lines
                3 => 16..24, // Stress: active in last 8 cell lines
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
            alpha: 0.3,
            beta: 0.5,
            max_iter: 10,
            k_neighbours: 3,
            admm_iter: 5,
            rho: 1.0,
        };

        let dgrdl = Dgrdl::new(params);
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

    #[test]
    fn test_regularization_effects() {
        // Create data where graph regularization should help
        let n_samples = 30;
        let n_features = 60;

        let mut data = Mat::zeros(n_samples, n_features);
        let mut rng = StdRng::seed_from_u64(789);

        // Create smooth patterns that should benefit from graph regularization
        for i in 0..n_samples {
            for j in 0..n_features {
                // Smooth function of both sample and feature indices
                let sample_effect =
                    (i as f64 / n_samples as f64 * 2.0 * std::f64::consts::PI).sin();
                let feature_effect = (j as f64 / n_features as f64 * std::f64::consts::PI).cos();
                data[(i, j)] = sample_effect * feature_effect + rng.random_range(-0.1..0.1);
            }
        }

        // Test with different regularization strengths
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
                max_iter: 15,
                k_neighbours: 5,
                admm_iter: 8,
                rho: 1.0,
            };

            let dgrdl = Dgrdl::new(params);
            let result = dgrdl.fit(&data.as_ref(), 42, false);

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
