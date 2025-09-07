use extendr_api::prelude::*;

use crate::core::base::pca_svd::*;
use crate::core::base::utils::scale_matrix_col;
use crate::utils::general::nested_vector_to_faer_mat;
use crate::utils::r_rust_interface::*;

/// Rust implementation of prcomp
///
/// @description Runs the singular value decomposition over the matrix x.
/// Assumes that samples = rows, and columns = features.
///
/// @param x Numeric matrix. Rows = samples, columns = features.
/// @param scale Boolean. Shall the columns additionally be scaled.
///
/// @return A list with:
/// \itemize{
///   \item scores - The product of x (centred and potentially scaled) with v.
///   \item v - v matrix of the SVD.
///   \item s - Eigenvalues of the SVD.
///   \item scaled - Boolean. Was the matrix scaled.
/// }
///
/// @export
#[extendr]
fn rs_prcomp(x: RMatrix<f64>, scale: bool) -> List {
    let x = r_matrix_to_faer(&x);
    let x_scaled = scale_matrix_col(&x.as_ref(), scale);
    let nrow = x_scaled.nrows() as f64;
    let svd_res = x_scaled.thin_svd().unwrap();
    let scores = x_scaled * svd_res.V();
    let s = svd_res
        .S()
        .column_vector()
        .iter()
        .cloned()
        .collect::<Vec<f64>>();
    let sdev: Vec<f64> = s.iter().map(|x| x / (nrow as f64 - 1.0).sqrt()).collect();

    list!(
        scores = faer_to_r_matrix(scores.as_ref()),
        v = faer_to_r_matrix(svd_res.V().as_ref()),
        s = sdev,
        scaled = scale,
    )
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
    let random_svd_res = randomised_svd(x, rank, seed, oversampling, n_power_iter);

    list!(
        u = faer_to_r_matrix(random_svd_res.u.as_ref()),
        v = faer_to_r_matrix(random_svd_res.v.as_ref()),
        s = random_svd_res.s
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
    return_loadings: bool,
) -> List {
    let target_covar = r_matrix_to_faer(&target_covar);
    let background_covar = r_matrix_to_faer(&background_covar);
    let target_mat = r_matrix_to_faer(&target_mat);

    let final_covar = target_covar - alpha * background_covar;

    let cpca_results = get_top_eigenvalues(&final_covar, n_pcs);

    let eigenvectors: Vec<Vec<f64>> = cpca_results.iter().map(|x| x.1.clone()).collect();

    let c_pca_loadings = nested_vector_to_faer_mat(eigenvectors, true);

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

extendr_module! {
    mod r_svd_pca;
    fn rs_prcomp;
    fn rs_random_svd;
    fn rs_contrastive_pca;
}
