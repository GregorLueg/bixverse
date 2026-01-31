use rayon::iter::*;

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

////////////
// Others //
////////////
