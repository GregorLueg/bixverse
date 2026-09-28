use bixverse_rs::core::math::matrix_helpers::scale_matrix_col;
use bixverse_rs::core::math::pca_missing::*;
use bixverse_rs::core::math::pca_svd::*;
use bixverse_rs::prelude::*;
use bixverse_rs::utils::matrix_utils::nested_vector_to_faer_mat;
use extendr_api::prelude::*;
use faer::Mat;

/////////////
// extendR //
/////////////

extendr_module! {
    mod r_svd_pca;
    fn rs_prcomp;
    fn rs_random_svd;
    fn rs_contrastive_pca;
    fn rs_ppca;
    fn rs_bpca;
}

//////////
// SVDs //
//////////

/// Rust implementation of prcomp
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Runs the singular value decomposition over the matrix x. Assumes that
/// samples = rows, and columns = features.
///
/// @param x Numeric matrix. Rows = samples, columns = features.
/// @param scale Boolean. Shall the columns be variance normalised. (Mean
/// centring will automatically occur.)
/// @param top_pcs Optional integer. Only return the top PCs (under the hood
/// all of them will be calculated). `NULL` returns all.
///
/// @returns A list with:
/// \itemize{
///   \item scores - The product of x (centred and potentially scaled) with v.
///   \item v - v matrix of the SVD.
///   \item s - Standard deviations of the PCs, i.e. singular values divided
///   by `sqrt(nrow(x) - 1)`.
///   \item scaled - Boolean. Was the matrix scaled.
/// }
///
/// @export
#[extendr]
fn rs_prcomp(
    x: RMatrix<f64>,
    scale: bool,
    top_pcs: Option<usize>,
) -> Result<List, extendr_api::Error> {
    let x = r_matrix_to_faer(&x);
    let x_scaled = scale_matrix_col(&x.as_ref(), scale);
    let nrow = x_scaled.nrows() as f64;
    let svd_res = x_scaled
        .thin_svd()
        .map_err(|e| BixverseErrors::FaerSvdError(format!("{e:?}")))
        .to_extendr()?;
    let scores = x_scaled * svd_res.V();
    let n = top_pcs.unwrap_or(scores.ncols());
    let sdev: Vec<f64> = svd_res
        .S()
        .column_vector()
        .iter()
        .take(n)
        .map(|x| x / (nrow - 1.0).sqrt())
        .collect();

    Ok(list!(
        scores = faer_to_r_matrix(scores.as_ref().subcols(0, n)),
        v = faer_to_r_matrix(svd_res.V().subcols(0, n)),
        s = sdev,
        scaled = scale,
    ))
}

/// Run randomised SVD over a matrix
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Runs a randomised singular value decomposition over a matrix. This
/// implementation is faster than the full SVD on large data sets, with slight
/// loss in precision.
///
/// @param x Numeric matrix. Rows = samples, columns = features.
/// @param scale Boolean. Shall the columns be variance normalised. (Mean
/// centring will automatically occur.)
/// @param rank Integer. The rank to use.
/// @param seed Integer. Random seed for reproducibility.
/// @param oversampling Optional integer. Defaults to `10L` if `NULL`.
/// @param n_power_iter Optional integer. Number of power iterations (each with
/// a QR decomposition). Defaults to `2L` if `NULL`.
///
/// @returns A list with:
/// \itemize{
///   \item scores - u matrix of the SVD multiplied by the singular values.
///   \item v - v matrix of the SVD.
///   \item s - Singular values of the SVD.
///   \item scaled - Boolean. Was the matrix scaled.
/// }
///
/// @export
#[extendr]
fn rs_random_svd(
    x: RMatrix<f64>,
    scale: bool,
    rank: usize,
    seed: usize,
    oversampling: Option<usize>,
    n_power_iter: Option<usize>,
) -> extendr_api::Result<List> {
    let x = r_matrix_to_faer(&x);
    let x_scaled = scale_matrix_col(&x.as_ref(), scale);
    let res =
        randomised_svd(x_scaled.as_ref(), rank, seed, oversampling, n_power_iter).to_extendr()?;
    let scores = Mat::<f64>::from_fn(res.u.nrows(), rank, |i, j| res.u[(i, j)] * res.s[j]);

    Ok(list!(
        scores = faer_to_r_matrix(scores.as_ref()),
        v = faer_to_r_matrix(res.v.as_ref()),
        s = res.s,
        scaled = scale
    ))
}

//////////////////////
// Constrastive PCA //
//////////////////////

/// Calculate the contrastive PCA
///
/// @description
/// `r lifecycle::badge("experimental")`
/// This function calculates the contrastive PCA given a target covariance
/// matrix and the background covariance matrix you wish to subtract. The alpha
/// parameter controls how much of the background covariance you wish to remove.
/// You have the options to return the feature loadings and you can specify
/// the number of cPCAs to return.
///
/// @param target_covar The co-variance matrix of the target data set.
/// @param background_covar The co-variance matrix of the background data set.
/// @param target_mat The original values of the target matrix. Rows =
/// samples, columns = features.
/// @param alpha How much of the background co-variance should be removed.
/// @param n_pcs How many contrastive PCs to return
/// @param return_loadings Shall the loadings be returned from the contrastive
/// PCA
///
/// @returns A list containing:
///  \itemize{
///   \item factors - The factors of the contrastive PCA, i.e. `target_mat`
///    multiplied by the loadings. Samples x `n_pcs`.
///   \item loadings - The loadings (top eigenvectors) of the contrastive PCA.
///    Features x `n_pcs`. Will be `NULL` if `return_loadings = FALSE`.
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
) -> Result<List, extendr_api::Error> {
    let target_covar = r_matrix_to_faer(&target_covar);
    let background_covar = r_matrix_to_faer(&background_covar);
    let target_mat = r_matrix_to_faer(&target_mat);

    let final_covar = target_covar - alpha * background_covar;

    let cpca_results = get_top_eigenvalues(&final_covar, n_pcs).to_extendr()?;

    let eigenvectors: Vec<Vec<f64>> = cpca_results.iter().map(|x| x.1.clone()).collect();

    let c_pca_loadings = nested_vector_to_faer_mat(eigenvectors, true);

    let c_pca_factors = target_mat * c_pca_loadings.clone();

    if return_loadings {
        Ok(list!(
            factors = faer_to_r_matrix(c_pca_factors.as_ref()),
            loadings = faer_to_r_matrix(c_pca_loadings.as_ref())
        ))
    } else {
        Ok(list!(
            factors = faer_to_r_matrix(c_pca_factors.as_ref()),
            loadings = r!(NULL)
        ))
    }
}

/////////////////////////////
// PCA with missing values //
/////////////////////////////

/// Turn the missing-value PCA results into an R list.
fn missing_pca_to_r_list(res: MissingPcaResults<f64>) -> List {
    list!(
        scores = faer_to_r_matrix(res.scores.as_ref()),
        loadings = faer_to_r_matrix(res.loadings.as_ref()),
        r2_cum = res.r2_cum,
        centre = res.centre,
        scale = res.scale,
        completed = faer_to_r_matrix(res.completed.as_ref()),
        noise_var = res.noise_var,
        n_iter = res.n_iter as i32,
        converged = res.converged
    )
}

/// Probabilistic PCA on a matrix with missing values
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Port of `ppca()` from pcaMethods. Fits the principal subspace by EM on the
/// observed entries only and imputes the missing ones from it. The start is
/// drawn from `seed` with the Rust RNG, so results match pcaMethods at
/// convergence, not iterate by iterate.
///
/// @param x Numeric matrix. Rows = samples, columns = features. `NA` marks a
/// missing value.
/// @param ppca_params List. The PPCA parameters, see
/// [bixverse::params_ppca()].
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @returns A list with:
/// \itemize{
///   \item scores - Samples x `n_pcs`.
///   \item loadings - Features x `n_pcs`, orthonormal.
///   \item r2_cum - Cumulative R^2 per component on the completed matrix.
///   \item centre - Column centres that were subtracted.
///   \item scale - Column scales that were divided out.
///   \item completed - `x` with the missing entries imputed.
///   \item noise_var - Residual variance outside the subspace.
///   \item n_iter - EM iterations run.
///   \item converged - Whether `tol` was reached before `max_iter`.
/// }
///
/// @references Roweis, NIPS, 1998; Tipping and Bishop, J R Stat Soc B, 1999;
/// Stacklies, et al., Bioinformatics, 2007
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_ppca(x: RMatrix<f64>, ppca_params: List, verbose: usize) -> extendr_api::Result<List> {
    // R's NA_real_ is a NaN bit pattern, which is what the crate tests for
    // with is_nan(), so the matrix crosses without a copy
    let x = r_matrix_to_faer(&x);
    let params = PpcaParams::from_r_list(ppca_params)?;
    let res = ppca(x, Some(params), parse_verbosity_level(verbose)).to_extendr()?;

    Ok(missing_pca_to_r_list(res))
}

/// Bayesian PCA on a matrix with missing values
///
/// @description
/// `r lifecycle::badge("experimental")`
/// Port of `bpca()` from pcaMethods. Variational Bayes with an ARD prior per
/// component, so superfluous components shrink towards zero. Deterministic,
/// the start comes from an SVD.
///
/// @param x Numeric matrix. Rows = samples, columns = features. `NA` marks a
/// missing value.
/// @param bpca_params List. The BPCA parameters, see
/// [bixverse::params_bpca()].
/// @param verbose Integer. `0L` - quiet; `1L` - normal verbosity; `2L` -
/// detailed verbosity.
///
/// @returns A list with:
/// \itemize{
///   \item scores - Samples x `n_pcs`.
///   \item loadings - Features x `n_pcs`, not orthonormal.
///   \item r2_cum - Cumulative R^2 per component on the observed entries.
///   \item centre - Column centres that were subtracted.
///   \item scale - Column scales that were divided out.
///   \item completed - `x` with the missing entries imputed.
///   \item noise_var - Residual variance, `1 / tau`.
///   \item n_iter - Variational steps run.
///   \item converged - Whether `tol` was reached before `max_iter`.
/// }
///
/// @references Oba, et al., Bioinformatics, 2003; Stacklies, et al.,
/// Bioinformatics, 2007
///
/// @export
///
/// @keywords internal
#[extendr]
fn rs_bpca(x: RMatrix<f64>, bpca_params: List, verbose: usize) -> extendr_api::Result<List> {
    let x = r_matrix_to_faer(&x);
    let params = BpcaParams::from_r_list(bpca_params)?;
    let res = bpca(x, Some(params), parse_verbosity_level(verbose)).to_extendr()?;

    Ok(missing_pca_to_r_list(res))
}
