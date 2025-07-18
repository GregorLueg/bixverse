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
/// * `beta` - Gene effect regularisation weight
/// * `max_iter` - Maximum number of iterations
/// * `k_neighbors` - Number of neighbours in the KNN graph
/// * `admm_iter` - ADMM iterations for sparse coding
/// * `rho` - ADMM step size
#[derive(Debug, Clone)]
pub struct DgrdlParams {
    pub sparsity: usize,
    pub dict_size: usize,
    pub alpha: f64,
    pub beta: f64,
    pub max_iter: usize,
    pub k_neighbors: usize,
    pub admm_iter: usize,
    pub rho: f64,
}

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
    /// The `DgrdlParams` structure
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
        let k_neighbors = dgrdl_params
            .get("k_neighbors")
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
            k_neighbors,
            admm_iter,
            rho,
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
#[derive(Debug, Clone)]
pub struct Dgrdl {
    params: DgrdlParams,
}

impl Dgrdl {
    pub fn new(params: DgrdlParams) -> Self {
        Self { params }
    }

    // Run DGRDL algorithm
    pub fn fit(&self, data: &MatRef<f64>, verbose: bool) -> DgrdlResults {
        let n_features = data.ncols();

        let mut dictionary = self.initialise_dictionary(data);

        let feature_laplacian = get_dgrdl_laplacian(&data.as_ref(), self.params.k_neighbors, true);
        let sample_laplacian = get_dgrdl_laplacian(&data.as_ref(), self.params.k_neighbors, false);

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
                &feature_laplacian,
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

    fn initialise_dictionary(&self, data: &MatRef<f64>) -> Mat<f64> {
        let k = self.params.dict_size;
        let (n_samples, n_features) = data.shape();

        let mut medoids = Vec::new();
        let mut remaining: Vec<usize> = (0..n_features).collect();

        medoids.push(remaining.remove(0));

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
/// The Laplacian matrix
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
/// Alternating Direction Method of Multipliers (ADMM) algorithm to solve the
/// system
///
/// ### Params
///
/// ### Returns
pub fn grsc_admm(
    dictionary: &Mat<f64>,
    data: &MatRef<f64>,
    gene_laplacian: &Mat<f64>,
    sparsity: usize,
    beta: f64,
    max_iter: usize,
    rho: f64,
) -> Mat<f64> {
    let (_, k) = dictionary.shape();
    let (_, m) = data.shape();

    let gram = dictionary.transpose() * dictionary;
    let dt_y = dictionary.transpose() * data;

    // Initialize variables
    #[allow(unused_assignments)]
    let mut x: Mat<f64> = Mat::zeros(k, m);
    let mut z: Mat<f64> = Mat::zeros(k, m);
    let mut u: Mat<f64> = Mat::zeros(k, m);

    for _ in 0..max_iter {
        // X update: solve (G + ρI)X + βXL = D^T Y + ρ(Z - U)
        x = solve_x_update(&gram, &dt_y, &z, &u, gene_laplacian, beta, rho);

        // Z update: sparse projection
        z = sparse_projection(&(&x + &u), sparsity);

        // U update: dual variable
        u = &u + &x - &z;
    }

    z
}

fn solve_x_update(
    gram: &Mat<f64>,
    dt_y: &Mat<f64>,
    z: &Mat<f64>,
    u: &Mat<f64>,
    gene_laplacian: &Mat<f64>,
    beta: f64,
    rho: f64,
) -> Mat<f64> {
    let k = gram.nrows();

    let identity_mat: Mat<f64> = Mat::identity(k, k);
    let a: Mat<f64> = gram + &(rho * identity_mat);
    let b = beta * gene_laplacian;
    let c = dt_y + &(rho * (z - u));

    sylvester_solver(&a.as_ref(), &b.as_ref(), &c.as_ref())
}

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

fn column_distance(col_i: ColRef<f64>, col_j: ColRef<f64>) -> f64 {
    col_i
        .iter()
        .zip(col_j.iter())
        .map(|(a, b)| (a - b).powi(2))
        .sum::<f64>()
        .sqrt()
}

fn update_single_atom(
    data: &MatRef<f64>,
    coefficients: &Mat<f64>,
    atom_idx: usize,
    context_laplacian: &Mat<f64>,
    _gene_laplacian: &Mat<f64>,
    alpha: f64,
    _beta: f64,
) -> Mat<f64> {
    let n_contexts = data.nrows();

    // Find signals that use this atom
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

    // Simplified atom update (would implement full regularized version)
    let x_j_norm_sq = x_j_active.iter().map(|x| x * x).sum::<f64>();

    if x_j_norm_sq > 1e-12 {
        let regularization_term = alpha * context_laplacian;
        let identity: Mat<f64> = Mat::identity(n_contexts, n_contexts);
        let system_matrix = &(x_j_norm_sq * identity) + regularization_term;

        let rhs = residual * Mat::from_fn(x_j_active.len(), 1, |i, _| x_j_active[i]);

        let lu = PartialPivLu::new(system_matrix.as_ref());

        lu.solve(&rhs)
    } else {
        Mat::zeros(n_contexts, 1)
    }
}

fn update_dictionary(
    data: &MatRef<f64>,
    coefficients: &Mat<f64>,
    feature_laplacian: &Mat<f64>,
    sample_laplacian: &Mat<f64>,
    alpha: f64,
    beta: f64,
) -> Mat<f64> {
    let (n_contexts, k) = data.shape();
    let mut dictionary = Mat::zeros(n_contexts, k);

    // Update each dictionary atom
    for atom_idx in 0..k {
        let updated_atom = update_single_atom(
            data,
            coefficients,
            atom_idx,
            feature_laplacian,
            sample_laplacian,
            alpha,
            beta,
        );

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
