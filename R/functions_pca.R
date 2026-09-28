# pca with missing values ------------------------------------------------------

## helpers ---------------------------------------------------------------------

#' Carry the input dimnames over to a missing-value PCA result
#'
#' @param res List. Output of [bixverse::rs_ppca()] or
#' [bixverse::rs_bpca()].
#' @param x Numeric matrix. The input matrix.
#'
#' @returns `res` with row names on `scores`, row names on `loadings` and the
#' dimnames of `x` on `completed`. Components are named `PC1`, `PC2`, etc.
#'
#' @keywords internal
.name_missing_pca_res <- function(res, x) {
  checkmate::assertList(res)
  checkmate::assertMatrix(x)

  pc_names <- sprintf("PC%i", seq_len(ncol(res$scores)))
  dimnames(res$scores) <- list(rownames(x), pc_names)
  dimnames(res$loadings) <- list(colnames(x), pc_names)
  dimnames(res$completed) <- dimnames(x)

  res
}

#' Assert that no row or column is entirely missing
#'
#' @description
#' Checked on the R side so the error names 1-based indices.
#'
#' @param x Numeric matrix.
#'
#' @returns Invisibly `x`, or an error naming the empty rows or columns.
#'
#' @keywords internal
.assert_observed_margins <- function(x) {
  checkmate::assertMatrix(x)

  observed <- !is.na(x)
  empty_rows <- which(rowSums(observed) == 0L)
  empty_cols <- which(colSums(observed) == 0L)
  if (length(empty_rows) > 0L) {
    stop(sprintf(
      "Rows entirely missing: %s. Remove them first.",
      paste(empty_rows, collapse = ", ")
    ))
  }
  if (length(empty_cols) > 0L) {
    stop(sprintf(
      "Columns entirely missing: %s. Remove them first.",
      paste(empty_cols, collapse = ", ")
    ))
  }

  invisible(x)
}

## main functions --------------------------------------------------------------

#' Probabilistic PCA on a matrix with missing values
#'
#' @description
#' `r lifecycle::badge("experimental")`
#' Port of `ppca()` from pcaMethods. Fits the principal subspace by EM on the
#' observed entries only, then imputes the missing entries from it. Use it to
#' get a PCA or a completed matrix out of data with sporadic `NA`s, for
#' example proteomics or metabolomics intensities.
#'
#' @details
#' The initial loadings are drawn from `ppca_params$seed` with the Rust RNG, so
#' results agree with pcaMethods at convergence, not iterate by iterate. PPCA
#' loadings are orthonormal. Signs of the components are arbitrary.
#'
#' @param x Numeric matrix. Rows = samples, columns = features. `NA` (or
#' `NaN`) marks a missing value. No row or column may be entirely missing and
#' observed values must be finite.
#' @param ppca_params List, see [bixverse::params_ppca()]. `n_pcs` must not
#' exceed `min(dim(x))`.
#' @param .verbose Boolean or integer. Verbosity.
#'
#' @returns A list with:
#' \itemize{
#'   \item scores - Samples x `n_pcs`.
#'   \item loadings - Features x `n_pcs`, orthonormal.
#'   \item r2_cum - Cumulative R^2 per component on the completed matrix.
#'   \item centre - Column centres that were subtracted.
#'   \item scale - Column scales that were divided out.
#'   \item completed - `x` with the missing entries imputed. Observed entries
#'   are untouched.
#'   \item noise_var - Residual variance outside the subspace.
#'   \item n_iter - EM iterations run.
#'   \item converged - Whether `tol` was reached before `max_iter`.
#' }
#'
#' @references Roweis, NIPS, 1998; Tipping and Bishop, J R Stat Soc B, 1999;
#' Stacklies, et al., Bioinformatics, 2007
#'
#' @export
#'
#' @examples
#' set.seed(42L)
#' x <- matrix(rnorm(50L * 10L), nrow = 50L, ncol = 10L)
#' x[sample(length(x), 50L)] <- NA
#' ppca_res <- run_ppca(x, .verbose = FALSE)
#' head(ppca_res$scores)
run_ppca <- function(x, ppca_params = params_ppca(), .verbose = TRUE) {
  # checks
  checkmate::assertMatrix(x, mode = "numeric", min.rows = 2L, min.cols = 1L)
  checkmate::assertNumeric(x, finite = TRUE)
  .assert_observed_margins(x)
  assertPpcaParams(ppca_params)
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  storage.mode(x) <- "double"

  res <- rs_ppca(
    x = x,
    ppca_params = ppca_params,
    verbose = parse_verbosity(.verbose)
  )

  .name_missing_pca_res(res, x)
}

#' Bayesian PCA on a matrix with missing values
#'
#' @description
#' `r lifecycle::badge("experimental")`
#' Port of `bpca()` from pcaMethods. Variational Bayes with an ARD prior per
#' component, so components that the data do not support shrink towards zero
#' rather than fitting noise. Missing entries are imputed from the fit.
#'
#' @details
#' Deterministic, the start comes from an SVD. BPCA loadings are not
#' orthonormal, and `r2_cum` is computed on the observed entries only, as in
#' pcaMethods. Signs of the components are arbitrary.
#'
#' @inheritParams run_ppca
#'
#' @param bpca_params List, see [bixverse::params_bpca()]. `n_pcs` must not
#' exceed `min(dim(x))`.
#'
#' @returns A list with:
#' \itemize{
#'   \item scores - Samples x `n_pcs`.
#'   \item loadings - Features x `n_pcs`, not orthonormal.
#'   \item r2_cum - Cumulative R^2 per component on the observed entries.
#'   \item centre - Column centres that were subtracted.
#'   \item scale - Column scales that were divided out.
#'   \item completed - `x` with the missing entries imputed. Observed entries
#'   are untouched.
#'   \item noise_var - Residual variance, `1 / tau`.
#'   \item n_iter - Variational steps run.
#'   \item converged - Whether `tol` was reached before `max_iter`.
#' }
#'
#' @references Oba, et al., Bioinformatics, 2003; Stacklies, et al.,
#' Bioinformatics, 2007
#'
#' @export
#'
#' @examples
#' set.seed(42L)
#' x <- matrix(rnorm(50L * 10L), nrow = 50L, ncol = 10L)
#' x[sample(length(x), 50L)] <- NA
#' bpca_res <- run_bpca(x, .verbose = FALSE)
#' head(bpca_res$scores)
run_bpca <- function(x, bpca_params = params_bpca(), .verbose = TRUE) {
  # checks
  checkmate::assertMatrix(x, mode = "numeric", min.rows = 2L, min.cols = 1L)
  checkmate::assertNumeric(x, finite = TRUE)
  .assert_observed_margins(x)
  assertBpcaParams(bpca_params)
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  storage.mode(x) <- "double"

  res <- rs_bpca(
    x = x,
    bpca_params = bpca_params,
    verbose = parse_verbosity(.verbose)
  )

  .name_missing_pca_res(res, x)
}
