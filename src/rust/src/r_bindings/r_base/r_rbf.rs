use extendr_api::prelude::*;
use faer::Mat;

use crate::core::base::rbf::*;
use crate::utils::r_rust_interface::*;

/// Apply a Radial Basis Function
///
/// @description Applies a radial basis function (RBF) to a given distance
/// vector. Has at the option to apply a Gaussian, Bump or Inverse Quadratic
/// RBF.
///
/// @param x Numeric vector. The distances you wish to apply the Gaussian kernel
/// onto.
/// @param epsilon Float. Epsilon parameter for the RBF.
/// @param rbf_type String. Needs to be from `c("gaussian", "bump", "inverse_quadratic")`.
///
/// @return The affinities after the Kernel was applied.
///
/// @export
#[extendr]
fn rs_rbf_function(x: &[f64], epsilon: f64, rbf_type: &str) -> extendr_api::Result<Vec<f64>> {
    let rbf_fun = parse_rbf_types(rbf_type)
        .ok_or_else(|| extendr_api::Error::Other(format!("Invalid RBF function: {}", rbf_type)))?;

    let res: Vec<f64> = match rbf_fun {
        RbfType::Gaussian => rbf_gaussian(x, &epsilon),
        RbfType::Bump => rbf_bump(x, &epsilon),
        RbfType::InverseQuadratic => rbf_inverse_quadratic(x, &epsilon),
    };

    Ok(res)
}

/// Apply a Radial Basis Function (to a matrix)
///
/// @description Applies a radial basis function (RBF) to a given distance
/// matrix. Has at the option to apply a Gaussian, Bump or Inverse Quadratic
/// RBF.
///
/// @param x Numeric Matrix. The distances you wish to apply the Gaussian kernel
/// onto.
/// @param epsilon Float. Epsilon parameter for the RBF.
/// @param rbf_type String. Needs to be from `c("gaussian", "bump", "inverse_quadratic")`.
///
/// @return The affinities after the Kernel was applied.
///
/// @export
#[extendr]
fn rs_rbf_function_mat(
    x: RMatrix<f64>,
    epsilon: f64,
    rbf_type: &str,
) -> extendr_api::Result<extendr_api::RArray<f64, [usize; 2]>> {
    let x = r_matrix_to_faer(&x);

    let rbf_fun = parse_rbf_types(rbf_type)
        .ok_or_else(|| extendr_api::Error::Other(format!("Invalid RBF function: {}", rbf_type)))?;

    let res: Mat<f64> = match rbf_fun {
        RbfType::Gaussian => rbf_gaussian_mat(x, &epsilon),
        RbfType::Bump => rbf_bump_mat(x, &epsilon),
        RbfType::InverseQuadratic => rbf_inverse_quadratic_mat(x, &epsilon),
    };

    let res = faer_to_r_matrix(res.as_ref());

    Ok(res)
}

/// Helper to identify the right epsilon parameter
///
/// @description This function will take a distance vector from the upper
/// triangle of a symmetric distance matrix and apply the desired RBF with the
/// supplied epsilon from epsilon vec. Subsequently, the column sums will be
/// measured to identify the total similarity of each feature with other
/// features. This data can be used to see if the data follows scale-free
/// topology for example to identify the right epsilon parameter with the given
/// RBF.
///
/// @param dist Numeric vector. The distances you wish to apply the RBF function
/// to.
/// @param epsilon_vec Numeric vector. The epsilons you wish to use/test.
/// @param original_dim Integer. The original dimensions of the symmetric
/// distance matrix.
/// @param shift Integer. Was the matrix shifted up (0 = diagonal included; 1
/// diagonal not incldued).
/// @param rbf_type String. Option of `c('gaussian', 'bump')` for the currently
/// implemented RBF function.
///
/// @return A matrix with rows being the epsilons tested, and columns
/// representing the summed affinity to other features.
///
/// @export
#[extendr]
fn rs_rbf_iterate_epsilons(
    dist: &[f64],
    epsilon_vec: &[f64],
    original_dim: usize,
    shift: usize,
    rbf_type: &str,
) -> extendr_api::Result<extendr_api::RArray<f64, [usize; 2]>> {
    let band_width_data = rbf_iterate_epsilons(dist, epsilon_vec, original_dim, shift, rbf_type)?;

    Ok(faer_to_r_matrix(band_width_data.as_ref()))
}

extendr_module! {
    mod r_rbf;
    fn rs_rbf_function;
    fn rs_rbf_function_mat;
    fn rs_rbf_iterate_epsilons;
}
