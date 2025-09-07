use faer::{Mat, MatRef};
use rayon::iter::*;

use crate::core::base::utils::col_sums;
use crate::utils::general::{nested_vector_to_faer_mat, upper_triangle_to_sym_faer};

///////////
// Enums //
///////////

/// Enum for the RBF function
#[derive(Debug)]
pub enum RbfType {
    Gaussian,
    Bump,
    InverseQuadratic,
}

////////////
// Params //
////////////

/// Parsing the RBF function
///
/// ### Params
///
/// * `s` - String to transform into `RbfType`
///
/// ### Returns
///
/// Returns the `RbfType`
pub fn parse_rbf_types(s: &str) -> Option<RbfType> {
    match s.to_lowercase().as_str() {
        "gaussian" => Some(RbfType::Gaussian),
        "bump" => Some(RbfType::Bump),
        "inverse_quadratic" => Some(RbfType::InverseQuadratic),
        _ => None,
    }
}

///////////////////
// RBF functions //
///////////////////

/// Gaussian Radial Basis function
///
/// Applies a Gaussian Radial Basis function on a vector of distances with the
/// following formula:
/// ```
/// φ(r) = e^(-(εr)²)
/// ```
///
/// ### Params
///
/// * `dist` - Vector of distances
/// * `epsilon` - Shape parameter controlling function width
///
/// ### Returns
///
/// The resulting affinity vector
pub fn rbf_gaussian(dist: &[f64], epsilon: &f64) -> Vec<f64> {
    dist.par_iter()
        .map(|x| f64::exp(-((x * *epsilon).powi(2))))
        .collect()
}

/// Gaussian Radial Basis function for matrices.
///
/// Applies a Gaussian Radial Basis function on a matrix of distances with the
/// following formula:
///
/// ```
/// φ(r) = e^(-(εr)²)
/// ```
///
/// ### Params
///
/// * `dist` - Matrix of distances
/// * `epsilon` - Shape parameter controlling function width
///
/// ### Returns
///
/// The affinity matrix
pub fn rbf_gaussian_mat(dist: MatRef<f64>, epsilon: &f64) -> Mat<f64> {
    let ncol = dist.ncols();
    let nrow = dist.nrows();
    Mat::from_fn(nrow, ncol, |i, j| {
        let x = dist.get(i, j);
        f64::exp(-((x * epsilon).powi(2)))
    })
}

/// Bump Radial Basis function
///
/// Applies a Bump Radial Basis function on a vector of distances with the
/// following formula:
/// ```
/// φ(r) = { exp(-1/(1-(εr)²)) + 1,  if εr < 1
///        { 0,                      if εr ≥ 1
/// ```
///
/// ### Params
///
/// * `dist` - Vector of distances
/// * `epsilon` - Shape parameter controlling function width
///
/// ### Returns
///
/// The resulting affinity vector
pub fn rbf_bump(dist: &[f64], epsilon: &f64) -> Vec<f64> {
    dist.par_iter()
        .map(|x| {
            if *x < (1.0 / epsilon) {
                f64::exp(-(1_f64 / (1_f64 - (*epsilon * x).powi(2))) + 1_f64)
            } else {
                0_f64
            }
        })
        .collect()
}

/// Bump Radial Basis function for matrices
///
/// Applies a Bump Radial Basis function on a matrix of distances with the
/// following formula:
/// ```
/// φ(r) = { exp(-1/(1-(εr)²)) + 1,  if εr < 1
///        { 0,                      if εr ≥ 1
/// ```
///
/// ### Params
///
/// * `dist` - Matrix of distances
/// * `epsilon` - Shape parameter controlling function width
///
/// ### Returns
///
/// The resulting affinity matrix
pub fn rbf_bump_mat(dist: MatRef<f64>, epsilon: &f64) -> Mat<f64> {
    let ncol = dist.ncols();
    let nrow = dist.nrows();
    Mat::from_fn(nrow, ncol, |i, j| {
        let x = dist.get(i, j);
        if *x < (1.0 / epsilon) {
            f64::exp(-(1_f64 / (1_f64 - (*epsilon * x).powi(2))) + 1_f64)
        } else {
            0_f64
        }
    })
}

/// Inverse quadratic RBF
///
/// Applies a Inverse Quadratic Radial Basis function on a vector of distances
/// with the following formula:
/// ```
/// φ(r) = 1/(1 + (εr)²)
/// ```
///
/// ### Params
///
/// * `dist` - Vector of distances
/// * `epsilon` - Shape parameter controlling function width
///
/// ### Return
///
/// The resulting affinity vector
pub fn rbf_inverse_quadratic(dist: &[f64], epsilon: &f64) -> Vec<f64> {
    dist.par_iter()
        .map(|x| 1.0 / (1.0 + (*epsilon * x).powi(2)))
        .collect()
}

/// Inverse quadratic RBF for matrices
///
/// Applies a Inverse Quadratic Radial Basis function on a matrix of distances
/// with the following formula:
/// ```
/// φ(r) = 1/(1 + (εr)²)
/// ```
///
/// ### Params
///
/// * `dist` - Matrix of distances
/// * `epsilon` - Shape parameter controlling function width
///
/// ### Returns
///
/// The resulting affinity matrix
pub fn rbf_inverse_quadratic_mat(dist: MatRef<f64>, epsilon: &f64) -> Mat<f64> {
    let ncol = dist.ncols();
    let nrow = dist.nrows();
    Mat::from_fn(nrow, ncol, |i, j| {
        let x = dist.get(i, j);
        1.0 / (1.0 + (*epsilon * x).powi(2))
    })
}

////////////
// Others //
////////////

/// Test different epsilons over a distance vector
///
/// ### Params
///
/// * `dist` - The distance vector on which to apply the specified RBF function.
///   Assumes that these are the values of upper triangle of the distance
///   matrix.
/// * `epsilons` - Vector of epsilons to test.
/// * `n` - Original dimensions of the distance matrix from which `dist` was
///   derived.
/// * `shift` - Was a shift applied during the generation of the vector, i.e., was
///   the diagonal included or not.
/// * `rbf_type` - Which RBF function to apply on the distance vector.
///
/// ### Returns
///
/// The column sums of the resulting adjacency matrices after application of the
/// RBF function to for example check if these are following power law distributions.
pub fn rbf_iterate_epsilons(
    dist: &[f64],
    epsilons: &[f64],
    n: usize,
    shift: usize,
    rbf_type: &str,
) -> Result<faer::Mat<f64>, String> {
    // Now specifying String as the error type
    let rbf_fun =
        parse_rbf_types(rbf_type).ok_or_else(|| format!("Invalid RBF function: {}", rbf_type))?;

    let k_res: Vec<Vec<f64>> = epsilons
        .par_iter()
        .map(|epsilon| {
            let affinity_adj = match rbf_fun {
                RbfType::Gaussian => rbf_gaussian(dist, epsilon),
                RbfType::Bump => rbf_bump(dist, epsilon),
                RbfType::InverseQuadratic => rbf_inverse_quadratic(dist, epsilon),
            };
            let affinity_adj_mat = upper_triangle_to_sym_faer(&affinity_adj, shift, n);
            col_sums(affinity_adj_mat.as_ref())
        })
        .collect();

    Ok(nested_vector_to_faer_mat(k_res, true))
}
