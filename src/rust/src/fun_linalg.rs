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
/// @return An R list with loadings and factors. If return_loadings == FALSE, 
/// then loadings will be NULL.
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
      loadings = faer_to_r_matrix(c_pca_loadings),
      factors = faer_to_r_matrix(c_pca_factors)
    )
  } else {
    list!(
      loadings = r!(NULL),
      factors = faer_to_r_matrix(c_pca_factors)
    )
  }
}


/// @export
#[extendr]
fn rs_whiten_matrix(
  x: RMatrix<f64>,
) -> extendr_api::RArray<f64, [usize; 2]> {
  let x = r_matrix_to_faer(x);
  
  let whiten = whiten_matrix(x);

  faer_to_r_matrix(whiten)
}


extendr_module! {
  mod fun_linalg;
  fn rs_covariance;
  fn rs_contrastive_pca;
  fn rs_whiten_matrix;
}