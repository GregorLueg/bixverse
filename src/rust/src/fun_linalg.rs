use extendr_api::prelude::*;

use crate::helpers_linalg::*;
use crate::utils_r_rust::{r_matrix_to_faer, faer_to_r_matrix};
use crate::utils_rust::{nested_vector_to_faer_mat, upper_triangle_indices};


/// Calculate the column-wise co-variance.
/// 
/// @description Calculates the co-variance of the columns.
/// WARNING! Incorrect use can cause kernel crashes. Wrapper around the Rust 
/// functions with type checks are provided in the package.
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
  let mat = r_matrix_to_faer(&x);
  let covar = column_covariance(&mat);

  faer_to_r_matrix(covar.as_ref())
}

/// Calculate the column wise correlations.
/// 
/// @description Calculates the correlation matrix of the columns.
/// WARNING! Incorrect use can cause kernel crashes. Wrapper around the Rust 
/// functions with type checks are provided in the package.
/// 
/// @param x R matrix with doubles.
/// @param spearman Shall the Spearman correlation be calculated instead of 
/// Pearson.
/// 
/// @returns The correlation matrix.
/// 
/// @export
#[extendr]
fn rs_cor(
  x: RMatrix<f64>,
  spearman: bool
) -> extendr_api::RArray<f64, [usize; 2]> {
  let mat = r_matrix_to_faer(&x);

  let cor = column_correlation(
    mat,
    spearman
  );

  faer_to_r_matrix(cor.as_ref())
}

/// Calculate the column wise correlations.
/// 
/// @description Calculates the correlation matrix of the columns. This function
/// will return the upper triangle. WARNING! Incorrect use can cause kernel 
/// crashes. Wrapper around the Rust functions with type checks are provided in
/// the package.
/// 
/// @param x R matrix with doubles.
/// @param spearman Shall the Spearman correlation be calculated instead of 
/// Pearson.
/// @param shift Shall a shift be applied to the matrix. 0 = the diagonal will
/// be included. 1 = the diagonal will not be included.
/// 
/// @returns The upper triangle of the correlation matrix iterating through the
/// rows, shifted by one (the diagonal will not be returned).
/// 
/// @export
#[extendr]
fn rs_cor_upper_triangle(
  x: RMatrix<f64>,
  spearman: bool,
  shift: usize,
) -> Vec<f64> {
  // Calculate the correlations
  let mat = r_matrix_to_faer(&x);
  let cor = column_correlation(
    mat,
    spearman
  );
  let upper_triangle_indices = upper_triangle_indices(
    mat.ncols(),
    shift
  );
  let mut cor_flat = Vec::new();
  for (&r, &c) in upper_triangle_indices.0
    .iter()
    .zip(upper_triangle_indices.1.iter()) {
    cor_flat.push(*cor.get(r, c));
  }

  cor_flat
}

/// Calculate the column wise differential correlation between two sets of data.
/// 
/// @description This function calculates the differential correlation based on
/// the Fisher method. For speed purposes, the function will only calculate the
/// differential correlation on the upper triangle of the two correlation
/// matrices.
/// WARNING! Incorrect use can cause kernel crashes. Wrapper around the Rust 
/// functions with type checks are provided in the package.
/// 
/// @param x_a R matrix a to be used for the differential correlation analysis.
/// @param x_b R matrix a to be used for the differential correlation analysis.
/// @param spearman Shall the Spearman correlation be calculated instead of 
/// Pearson.
/// 
/// @return A list containing:
///  \itemize{
///   \item r_a - The correlation coefficients in the upper triangle of 
///   matrix a.
///   \item r_b - The correlation coefficients in the upper triangle of 
///   matrix b.
///   \item z_score - The z-scores of the difference in correlation 
///   coefficients. 
///   \item p_val - The z-scores transformed to p-values.
/// }
/// 
/// @export
#[extendr]
fn rs_differential_cor(
  x_a: RMatrix<f64>,
  x_b: RMatrix<f64>,
  spearman: bool
) -> List {
  assert!(
    x_a.ncols() == x_b.ncols(),
    "Input matrices must have the same number of columns. Found {} columns in first matrix and {} in second.",
    x_a.ncols(), 
    x_b.ncols()
  );
  let n_sample_a = x_a.nrows();
  let n_sample_b = x_b.nrows();
  let mat_a = r_matrix_to_faer(&x_a);
  let mat_b = r_matrix_to_faer(&x_b);

  let cor_a = column_correlation(mat_a, spearman);
  let cor_b = column_correlation(mat_b, spearman);

  let diff_cor = calculate_diff_correlation(
    &cor_a,
    &cor_b,
    n_sample_a,
    n_sample_b,
    spearman
  );

  list!(
    r_a = diff_cor.r_a,
    r_b = diff_cor.r_b,
    z_score = diff_cor.z_score,
    p_val = diff_cor.p_vals
  )
}


/// Calculate the contrastive PCA
/// 
/// @description This function calculate the contrastive PCA given a target
/// covariance matrix and the background covariance matrix you wish to subtract.
/// The alpha parameter controls how much of the background covariance you wish
/// to remove. You have the options to return the feature loadings and you can
/// specificy the number of cPCAs to return. WARNING! Incorrect use can cause
/// kernel crashes. Wrapper around the Rust functions with type checks are
/// provided in the package.
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
  let target_covar = r_matrix_to_faer(&target_covar);
  let background_covar = r_matrix_to_faer(&background_covar);
  let target_mat = r_matrix_to_faer(&target_mat);

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
      factors = faer_to_r_matrix(c_pca_factors.as_ref()),
      loadings = faer_to_r_matrix(c_pca_loadings.as_ref())
    )
  } else {
    list!(
      factors = faer_to_r_matrix(c_pca_factors.as_ref()),
      loadings = r!(NULL)
    )
  }
}


/// Run randomised SVD over a matrix
/// 
/// @description Runs a randomised singular value decomposition over a matrix.
/// This implementation is faster than the full SVD on large data sets, with 
/// slight loss in precision.
/// 
/// @param x Numeric matrix. Rows = samples, columns = features.
/// @param rank Integer. The rank to use.
/// @param seed Integer. Random seed for reproducibility.
/// @param oversampling Integer. Defaults to `10L` if nothing is provided.
/// @param n_power_iter Integer. How often shall the QR decomposition be 
/// applied. Defaults to `2L` if nothing is provided.
/// 
/// @return A list with:
/// \itemize{
///   \item u - u matrix of the SVD.
///   \item v - v matrix of the SVD.
///   \item s - Eigenvalues of the SVD.
/// }
/// 
/// @export
#[extendr]
fn rs_random_svd(
  x: RMatrix<f64>,
  rank: usize,
  seed: usize,
  oversampling: Option<usize>,
  n_power_iter: Option<usize>,
) -> List {
  let x = r_matrix_to_faer(&x);
  let random_svd_res = randomised_svd(
    x, 
    rank, 
    seed,
    oversampling,
    n_power_iter
  );

  list!(
    u = faer_to_r_matrix(random_svd_res.u.as_ref()),
    v = faer_to_r_matrix(random_svd_res.v.as_ref()),
    s = random_svd_res.s
  )
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
  let band_width_data = rbf_iterate_epsilons(
    dist,
    epsilon_vec,
    original_dim,
    shift,
    rbf_type
  )?;

  Ok(faer_to_r_matrix(band_width_data.as_ref()))
}


extendr_module! {
  mod fun_linalg;
  fn rs_covariance;
  fn rs_cor;
  fn rs_random_svd;
  fn rs_cor_upper_triangle;
  fn rs_rbf_iterate_epsilons;
  fn rs_differential_cor;
  fn rs_contrastive_pca;
}
