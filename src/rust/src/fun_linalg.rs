use extendr_api::prelude::*;

use crate::utils_r_rust::{r_matrix_to_faer, faer_to_r_matrix};
use crate::helpers_linalg::*;
use crate::utils_rust::nested_vector_to_faer_mat;


/// Calculate the column-wise co-variance
/// 
/// @param x R matrix with doubles.
/// 
/// @returns The co-variance matrix.
/// 
/// @export
#[extendr]
fn rs_covariance(
  x: RMatrix<f64>
) -> extendr_api::RArray<f64, [usize; 2]> {
  let mat = r_matrix_to_faer(x);
  let covar = column_covariance(&mat);

  faer_to_r_matrix(covar)
}



/// Calculate the contrastive PCA
/// 
/// @param target_covar The co-variance matrix of the target data set.
/// @param background_covar The co-variance matrix of the background data set.
/// @param target_mat The original values of the target matrix.
/// @param alpha How much of the background co-variance should be removed.
/// @param n_pcs How many contrastive PCs to return
/// @param return_loadings Shall the loadings be returned from the contrastive
/// PCA
/// 
/// @return A list containing:
///  \itemize{
///   \item factors - The factors of the contrastive PCA.
///   \item loadings - The loadings of the contrastive PCA. Will be NULL if 
///    return_loadings is set to FALSE.
/// }
/// 
/// @export
#[extendr]
fn rs_contrastive_pca(
  target_covar: RMatrix<f64>,
  background_covar: RMatrix<f64>,
  target_mat: RMatrix<f64>,
  alpha: f64,
  n_pcs: usize,
  return_loadings: bool
) -> List {
  let target_covar = r_matrix_to_faer(target_covar);
  let background_covar = r_matrix_to_faer(background_covar);
  let target_mat = r_matrix_to_faer(target_mat);

  let final_covar = target_covar - alpha * background_covar;

  let cpca_results = get_top_eigenvalues(
    &final_covar, 
    n_pcs
  );

  let eigenvectors: Vec<Vec<f64>> = cpca_results
    .iter()
    .map(|x| { x.1.clone() })
    .collect();

  let c_pca_loadings= nested_vector_to_faer_mat(eigenvectors);

  let c_pca_factors = target_mat * c_pca_loadings.clone();

  if return_loadings {
    list!(
      factors = faer_to_r_matrix(c_pca_factors),
      loadings = faer_to_r_matrix(c_pca_loadings)
    )
  } else {
    list!(
      factors = faer_to_r_matrix(c_pca_factors),
      loadings = r!(NULL)
    )
  }
}

/// Whiten a matrix
/// 
/// @description Whitens the matrix for subsequent usage. This is a need pre-
/// processing step for ICA.
/// 
/// @param x The matrix to whiten. The whitening will happen over the columns.
/// 
/// @return A list containing:
///  \itemize{
///   \item whiten - The whitened matrix.
///   \item k - The K matrix.
/// }
/// 
/// @export
#[extendr]
fn rs_whiten_matrix(
  x: RMatrix<f64>,
) -> List {
  let x = r_matrix_to_faer(x);
  
  let (whiten, k) = whiten_matrix(x);

  list!(
    whiten = faer_to_r_matrix(whiten),
    k = faer_to_r_matrix(k)
  )
}

/// Run the Rust implementation of fast ICA.
/// 
/// @description This function serves as a wrapper over the fast ICA implementations
/// in Rust. It assumes a pre-whiten matrix and also an intialised w_init.
/// 
/// @param whiten The whitened matrix.
/// @param w_init The w_init matrix. ncols need to be equal to nrows of whiten.
/// @param maxit Maximum number of iterations to try if algorithm does not converge.
/// @param alpha The alpha parameter for the LogCosh implementation of ICA.
/// @param tol Tolerance parameter.
/// @param ica_type One of 'logcosh' or 'exp'.
/// @param verbose Controls the verbosity of the function.
/// 
/// @param x The matrix to whiten. The whitening will happen over the columns.
/// 
/// @return A list containing:
///  \itemize{
///   \item mixing - The mixing matrix for subsequent usage.
///   \item converged - Boolean if the algorithm converged.
/// }
/// 
/// @export
#[extendr]
fn rs_fast_ica(
  whiten: RMatrix<f64>,
  w_init: RMatrix<f64>,
  maxit: usize,
  alpha: f64,
  tol: f64,
  ica_type: &str,
  verbose: bool
) -> extendr_api::Result<List> {
  // assert!(!whiten.nrows() == w_init.ncols(), "The dimensions of the provided matrices don't work");

  let x = r_matrix_to_faer(whiten);
  let w_init = r_matrix_to_faer(w_init);

  let ica_type = parse_ica_type(ica_type).ok_or_else(|| format!("Invalid ICA type: {}", ica_type))?;

  let a = match ica_type {
    IcaType::Exp => fast_ica_exp(x, w_init, tol, maxit, verbose),
    IcaType::LogCosh => fast_ica_logcosh(x, w_init, tol, alpha, maxit, verbose),
  };

  Ok(list!(
    mixing = faer_to_r_matrix(a.0),
    converged = a.1 < tol)
  )
}

extendr_module! {
  mod fun_linalg;
  fn rs_covariance;
  fn rs_contrastive_pca;
  fn rs_whiten_matrix;
  fn rs_fast_ica;
}