# general ICA function ---------------------------------------------------------

#' Fast ICA via Rust
#'
#' @description
#' This functions is a wrapper over the Rust implementation of fastICA and the
#' generation of the pre-processed data and pre-whitening matrix. It has
#' the same two options `c("logcosh", "exp")` to run ICA in parallel modus. You
#' can control the parameters of ICA via `ica_params`.
#'
#' @param X Numeric matrix. The data on which you want to run fastICA.
#' @param n_icas Integer. Number of independent components to recover.
#' @param ica_fun String, element of `c("logcosh", "exp")`.
#' @param ica_params List. The ICA parameters, see
#' [bixverse::params_ica_general()] wrapper function. This function generates a
#' list containing:
#' \itemize{
#'  \item maxit - Integer. Maximum number of iterations for ICA.
#'  \item alpha - Float. The alpha parameter for the logcosh version of ICA.
#'  Should be between 1 to 2.
#'  \item max_tol - Maximum tolerance of the algorithm.
#'  \item verbose - Controls verbosity of the function.
#' }
#' @param fast_svd Boolean. Shall the randomised SVD be used. This is faster,
#' but less precise.
#' @param seed Integer. Seed to ensure reproducible results.
#'
#' @returns A list containing:
#' \itemize{
#'  \item w The mixing matrix w.
#'  \item A ICA results matrix A.
#'  \item S ICA results matrix S.
#'  \item converged Boolean indicating if algorithm converged.
#' }
#'
#' @export
#'
#' @importFrom zeallot %<-%
#'
#' @examples
#' # recover two mixed sources
#' sources <- cbind(sin((1:1000) / 20), rep(((1:200) - 100) / 100, 5))
#' mixed <- sources %*% matrix(c(0.291, 0.6557, -0.5439, 0.5572), 2, 2)
#' ica_res <- fast_ica_rust(mixed, n_icas = 2L, ica_fun = "logcosh", seed = 42L)
#' max(abs(cor(sources[, 1], t(ica_res$S))))
fast_ica_rust <- function(
  X,
  n_icas,
  ica_fun = c("logcosh", "exp"),
  ica_params = params_ica_general(),
  fast_svd = TRUE,
  seed = NULL
) {
  ica_fun <- match.arg(ica_fun)

  # scope
  X1 <- K <- NULL

  # checks
  checkmate::assertMatrix(X, mode = "numeric")
  checkmate::qassert(n_icas, sprintf("I1[2,%i]", ncol(X)))
  checkmate::assertChoice(ica_fun, c("logcosh", "exp"))
  checkmate::qassert(seed, c("I1", "0"))
  assertIcaParams(ica_params)

  c(X1, K) %<-%
    rs_prepare_whitening(
      x = X,
      fast_svd = fast_svd,
      seed = seed,
      rank = NULL,
      oversampling = NULL,
      n_power_iter = NULL
    )

  fast_ica_res <- fast_ica_rust_helper(
    X = X1,
    K = K,
    n_icas = n_icas,
    ica_fun = ica_fun,
    seed = seed
  )

  return(fast_ica_res)
}

## helpers ---------------------------------------------------------------------

#' Fast ICA via Rust from processed data
#'
#' @description
#' This functions is a wrapper over the Rust implementation of fastICA expected
#' already the pre-processed matrix and the pre-whitening matrix. It has
#' the same two options `c("logcosh", "exp")` to run ICA in parallel modus. You
#' can control the parameters of ICA via `ica_params`.
#'
#' @param X Numeric matrix. Processed data. Output of
#' [bixverse::rs_prepare_whitening()].
#' @param K The K matrix. Pre-whitening matrix. Output of
#' [bixverse::rs_prepare_whitening()].
#' @param n_icas Integer. Number of independent components to recover.
#' @param ica_fun String, element of `c("logcosh", "exp")`.
#' @param seed Integer. Seed to ensure reproducible results.
#' @param ica_params List. The ICA parameters, see [bixverse::params_ica_general()]
#' wrapper function. This function generates a list containing:
#' \itemize{
#'  \item maxit - Integer. Maximum number of iterations for ICA.
#'  \item alpha - Float. The alpha parameter for the logcosh version of ICA.
#'  Should be between 1 to 2.
#'  \item max_tol - Maximum tolerance of the algorithm.
#'  \item verbose - Controls verbosity of the function.
#' }
#'
#' @returns A list containing:
#' \itemize{
#'  \item w The mixing matrix w.
#'  \item A ICA results matrix A.
#'  \item S ICA results matrix S.
#'  \item converged Boolean indicating if algorithm converged.
#' }
#'
#' @export
#'
#' @importFrom zeallot %<-%
#'
#' @examples
#' # same run, but with the whitening done up front
#' sources <- cbind(sin((1:1000) / 20), rep(((1:200) - 100) / 100, 5))
#' mixed <- sources %*% matrix(c(0.291, 0.6557, -0.5439, 0.5572), 2, 2)
#' whitened <- rs_prepare_whitening(
#'   x = mixed,
#'   fast_svd = TRUE,
#'   seed = 42L,
#'   rank = NULL,
#'   oversampling = NULL,
#'   n_power_iter = NULL
#' )
#' ica_res <- fast_ica_rust_helper(
#'   X = whitened$x,
#'   K = whitened$k,
#'   n_icas = 2L,
#'   ica_fun = "logcosh",
#'   seed = 42L
#' )
#' dim(ica_res$S)
fast_ica_rust_helper <- function(
  X,
  K,
  n_icas,
  ica_fun = c("logcosh", "exp"),
  seed = NULL,
  ica_params = params_ica_general()
) {
  ica_fun <- match.arg(ica_fun)

  # Scope check
  a <- NULL
  # Checks
  checkmate::assertMatrix(X, mode = "numeric")
  checkmate::assertMatrix(K, mode = "numeric")
  checkmate::qassert(n_icas, sprintf("I1[2,%i]", ncol(X)))
  checkmate::assertChoice(ica_fun, c("logcosh", "exp"))
  checkmate::qassert(seed, c("I1", "0"))
  assertIcaParams(ica_params)
  # Function
  K <- matrix(K[1:n_icas, ], n_icas, dim(X)[1])
  X1 <- K %*% X
  set.seed(seed)
  w_init <- matrix(rnorm(n_icas * n_icas), nrow = n_icas, ncol = n_icas)

  c(a, converged) %<-%
    rs_fast_ica(X1, w_init, ica_type = ica_fun, ica_params = ica_params)

  w <- a %*% K
  S <- w %*% X
  A <- t(w) %*% solve(w %*% t(w))

  res <- list(
    w = w,
    A = A,
    S = S,
    converged = converged
  )

  return(res)
}
