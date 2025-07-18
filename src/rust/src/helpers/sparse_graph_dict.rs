use extendr_api::prelude::*;

use faer::{
    linalg::solvers::{PartialPivLu, Solve},
    ColRef, Mat, MatRef,
};

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
    pub fn fit(&self, data: &MatRef<f64>, verbose: bool) -> DgrdlResults {
        let n_features = data.ncols();

        let mut dictionary = self.initialise_dictionary(data);

        let feature_laplacian = get_dgrdl_laplacian(&data.as_ref(), self.params.k_neighbours, true);
        let sample_laplacian = get_dgrdl_laplacian(&data.as_ref(), self.params.k_neighbours, false);

        let mut coefficients: Mat<f64> = Mat::zeros(self.params.dict_size, n_features);

        for iter in 0..self.params.max_iter {
            if verbose {
                println!("DGRDL iteration {}/{}", iter + 1, self.params.max_iter);
            }

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
            dictionary = update_dictionary(
                data,
                &coefficients,
                &sample_laplacian,
                self.params.alpha,
                self.params.beta,
            );
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
    ///
    /// ### Returns
    ///
    /// Initialised dictionary matrix with normalized columns
    fn initialise_dictionary(&self, data: &MatRef<f64>) -> Mat<f64> {
        let k = self.params.dict_size;
        let (n_samples, n_features) = data.shape();

        let mut medoids = Vec::new();
        let mut remaining: Vec<usize> = (0..n_features).collect();

        // Start with first sample as initial medoid
        medoids.push(remaining.remove(0));

        // Greedy selection of remaining medoids
        for _ in 1..k {
            let mut best_medoid = 0;
            let mut best_cost = f64::INFINITY;

            for &candidate in &remaining {
                let cost = compute_medoid_cost(data, &medoids, candidate);
                if cost < best_cost {
                    best_cost = cost;
                    best_medoid = candidate;
                }
            }

            medoids.push(best_medoid);
            remaining.retain(|&x| x != best_medoid);
        }

        // Build dictionary from selected medoids
        let mut dictionary = Mat::zeros(n_samples, k);
        for (dict_idx, &gene_idx) in medoids.iter().enumerate() {
            let col = data.col(gene_idx);
            let norm = col.norm_l2();
            if norm > 1e-12 {
                for i in 0..n_samples {
                    dictionary[(i, dict_idx)] = col[i] / norm;
                }
            }
        }

        dictionary
    }
}

/////////////
// Helpers //
/////////////

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

    for j in 0..m {
        let col = x.col(j);
        let mut indexed_vals: Vec<(usize, f64)> =
            col.iter().enumerate().map(|(i, &val)| (i, val)).collect();

        indexed_vals.sort_by(|a, b| b.1.abs().partial_cmp(&a.1.abs()).unwrap());

        for &(idx, val) in indexed_vals.iter().take(sparsity) {
            result[(idx, j)] = val;
        }
    }

    result
}

/// Compute the cost of adding a candidate medoid
///
/// Calculates the total distance from all data points to their nearest
/// medoid (including the candidate), used for k-medoids clustering.
///
/// ### Params
///
/// * `data` - Input data matrix
/// * `current_medoids` - Currently selected medoid indices
/// * `candidate` - Candidate medoid index to evaluate
///
/// ### Returns
///
/// Total cost (sum of distances to nearest medoids) if candidate is added
fn compute_medoid_cost(data: &MatRef<f64>, current_medoids: &[usize], candidate: usize) -> f64 {
    let n_features = data.ncols();
    let mut total_cost = 0.0;

    for feature in 0..n_features {
        let mut min_dist = f64::INFINITY;

        // Distance to current medoids
        for &medoid in current_medoids {
            let dist = column_distance(data.col(feature), data.col(medoid));
            if dist < min_dist {
                min_dist = dist;
            }
        }

        // Distance to candidate
        let candidate_dist = column_distance(data.col(feature), data.col(candidate));
        if candidate_dist < min_dist {
            min_dist = candidate_dist;
        }

        total_cost += min_dist;
    }

    total_cost
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

/// Update a single dictionary atom
///
/// Updates one column of the dictionary by solving a regularized least squares
/// problem that incorporates sample graph regularization. Only considers
/// signals that actively use this atom (have non-zero coefficients).
///
/// ### Params
///
/// * `data` - Input data matrix Y
/// * `coefficients` - Current coefficient matrix X
/// * `atom_idx` - Index of the dictionary atom to update
/// * `context_laplacian` - Sample graph Laplacian L_s
/// * `alpha` - Sample regularization weight
/// * `_beta` - Feature regularization weight (unused in current implementation)
///
/// ### Returns
///
/// Updated dictionary atom as a column vector
fn update_single_atom(
    data: &MatRef<f64>,
    coefficients: &Mat<f64>,
    atom_idx: usize,
    sample_laplacian: &Mat<f64>,
    alpha: f64,
    _beta: f64,
) -> Mat<f64> {
    let n_contexts = data.nrows();

    // Find signals that use this atom (non-zero coefficients)
    let x_j = coefficients.row(atom_idx);
    let mut active_signals = Vec::new();
    let mut x_j_active = Vec::new();

    for (signal_idx, &coeff) in x_j.iter().enumerate() {
        if coeff.abs() > 1e-12 {
            active_signals.push(signal_idx);
            x_j_active.push(coeff);
        }
    }

    if active_signals.is_empty() {
        // Replace with random signal if unused
        return Mat::zeros(n_contexts, 1);
    }

    // Compute residual for active signals
    let mut residual = Mat::zeros(n_contexts, active_signals.len());
    for (res_idx, &signal_idx) in active_signals.iter().enumerate() {
        let signal = data.col(signal_idx);
        for i in 0..n_contexts {
            residual[(i, res_idx)] = signal[i];
        }
    }

    // Solve regularized least squares: (||x_j||²I + αL_s)d_j = residual * x_j
    let x_j_norm_sq = x_j_active.iter().map(|x| x * x).sum::<f64>();

    if x_j_norm_sq > 1e-12 {
        let regularization_term = alpha * sample_laplacian;
        let identity: Mat<f64> = Mat::identity(n_contexts, n_contexts);
        let system_matrix = &(x_j_norm_sq * identity) + regularization_term;

        let rhs = residual * Mat::from_fn(x_j_active.len(), 1, |i, _| x_j_active[i]);

        let lu = PartialPivLu::new(system_matrix.as_ref());

        lu.solve(&rhs)
    } else {
        Mat::zeros(n_contexts, 1)
    }
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
    beta: f64,
) -> Mat<f64> {
    let (n_contexts, _) = data.shape();
    let k = coefficients.nrows();
    let mut dictionary = Mat::zeros(n_contexts, k);

    // Update each dictionary atom
    for atom_idx in 0..k {
        let updated_atom =
            update_single_atom(data, coefficients, atom_idx, sample_laplacian, alpha, beta);

        // Normalize atom
        let norm = updated_atom.norm_l2();
        if norm > 1e-12 {
            for i in 0..n_contexts {
                dictionary[(i, atom_idx)] = updated_atom[(i, 0)] / norm;
            }
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
        let result = dgrdl.fit(&synthetic_data.as_ref(), true);

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
        let result = dgrdl.fit(&data.as_ref(), true);

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
        let result = dgrdl.fit(&data.as_ref(), true);

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
            let result = dgrdl.fit(&data.as_ref(), false);

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
