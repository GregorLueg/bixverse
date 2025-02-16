# helpers ----

#' Fast ICA via Rust
#'
#' @description
#' This functions is a wrapper over the Rust implementation of fastICA. It has
#' the same two options `c("logcosh", "exp")` to run ICA in the parallel modus.
#'
#' @param X_norm The scaled data. Output of [bixverse::rs_prepare_whitening()].
#' @param K The K matrix. Output of [bixverse::rs_prepare_whitening()].
#' @param n_icas Number of ICAs
#' @param ica_fun String, element of `c("logcosh", "exp")`.
#' @param seed Seed to ensure reproducible results.
#' @param maxit Maximum iterations to run. Defaults to 200L.
#' @param alpha Alpha parameters for the `"logcosh"` implementation. Defaults to
#' 1.
#' @param tol Tolerance at which the algorithm is considered converged. Defaults
#' to 1e-04.
#' @param .verbose Boolean controlling the verbosity.
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
#' @importFrom zeallot `%<-%`
fast_ica_rust <- function(X_norm,
                          K,
                          n_icas,
                          ica_fun = c("logcosh", "exp"),
                          seed = NULL,
                          maxit = 200L,
                          alpha = 1,
                          tol =  1e-04,
                          .verbose = FALSE) {
  # Checks
  checkmate::assertMatrix(X_norm, mode = "numeric")
  checkmate::assertMatrix(K, mode = "numeric")
  checkmate::qassert(n_icas, sprintf("I1[2,%i]", ncol(X_norm)))
  checkmate::assertChoice(ica_fun, c("logcosh", "exp"))
  checkmate::qassert(seed, c("I1", "0"))
  checkmate::qassert(.verbose, "B1")
  # Function
  K <- matrix(K[1:n_icas, ], n_icas, dim(X_norm)[1])
  X1 <- K %*% X_norm
  set.seed(seed)
  w_init <- matrix(rnorm(n_icas * n_icas), nrow = n_icas, ncol = n_icas)

  c(a, converged) %<-% rs_fast_ica(
    X1,
    w.init,
    maxit = maxit,
    alpha = alpha,
    tol = tol,
    ica_type = "logcosh",
    verbose = .verbose
  )

  w <- a %*% K
  S <- w %*% X_norm
  A <- t(w) %*% solve(w %*% t(w))

  res <- list(
    w <- w,
    A <- A,
    S <- S,
    converged <- converged
  )

  return(res)
}


# preprocessing ----

#' Prepare class for ICA
#'
#' @description
#' This is the generic function for doing the necessary preprocessing for
#' running independent component analysis.
#'
#' @export
ica_processing <- S7::new_generic(
  "ica_processing",
  "bulk_coexp"
)

#' @name ica_processing
#'
#' @description
#' This function will prepare the `bulk_coexp` for subsequent usage of the
#' contrastive ICA functions.
#'
#' @usage ...
#'
#' @param bulk_coexp `bulk_coexp` class, see [bixverse::bulk_coexp()].
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @return `bulk_coexp` with the needed data for ICA in the
#' properties of the class.
#'
#' @export
#'
#'
#' @importFrom magrittr `%>%`
#' @importFrom magrittr `%$%`
#' @importFrom zeallot `%<-%`
#' @import data.table
#'
#' @method ica_processing bulk_coexp
S7::method(ica_processing, bulk_coexp) <- function(bulk_coexp,
                                                   .verbose = TRUE) {
  # Checks
  checkmate::assertClass(bulk_coexp, "bixverse::bulk_coexp")
  checkmate::qassert(.verbose, "B1")

  # Function body
  if(purrr::is_empty(S7::prop(bulk_coexp, "processed_data")[['processed_data']])) {
    warning("No pre-processed data found. Defaulting to the raw data")
    target_mat <- S7::prop(bulk_coexp, "raw_data")
  } else {
    target_mat <- S7::prop(bulk_coexp, "processed_data")[['processed_data']]
  }

  # Transpose
  c(X1, K) %<-% rs_prepare_whitening(target_mat)

  S7::prop(bulk_coexp, "processed_data")[["X1"]] <- X1
  S7::prop(bulk_coexp, "processed_data")[["K"]] <- K

  return(bulk_coexp)
}
