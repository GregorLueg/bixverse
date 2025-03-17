# preprocessing ----------------------------------------------------------------

#' @title Prepare class for ICA
#'
#' @description
#' This is the generic function for doing the necessary preprocessing for
#' running independent component analysis.
#'
#' @param object The class, see [bixverse::bulk_coexp()]. Ideally, you
#' should run [bixverse::preprocess_bulk_coexp()] before applying this function.
#' @param fast_svd Boolean. Shall randomised SVD be used for the whitening.
#' This is faster and usually causes little precision loss.
#' @param random_seed Integer. Seed for the randomised SVD. Only relevant, if
#' fast_svd = `TRUE`.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @return `bulk_coexp` with the needed data for ICA in the
#' properties of the class.
#'
#' @export
ica_processing <- S7::new_generic(
  name = "ica_processing",
  dispatch_args = "object",
  fun = function(object,
                 fast_svd = TRUE,
                 random_seed = 123L,
                 .verbose = TRUE) {
    S7::S7_dispatch()
  }
)


#' @export
#'
#' @importFrom magrittr `%>%`
#' @importFrom magrittr `%$%`
#' @importFrom zeallot `%<-%`
#' @import data.table
#'
#' @method ica_processing bulk_coexp
S7::method(ica_processing, bulk_coexp) <- function(object,
                                                   fast_svd = TRUE,
                                                   random_seed = 123L,
                                                   .verbose = TRUE) {
  # Checks
  checkmate::assertClass(object, "bixverse::bulk_coexp")
  checkmate::qassert(.verbose, "B1")

  # Function body
  if(purrr::is_empty(S7::prop(object, "processed_data")[['processed_data']])) {
    warning("No pre-processed data found. Defaulting to the raw data")
    target_mat <- S7::prop(object, "raw_data")
  } else {
    target_mat <- S7::prop(object, "processed_data")[['processed_data']]
  }

  # Whiten the data
  if(.verbose)
    message("Preparing the whitening of the data for ICA.")
  if(fast_svd && .verbose)
    message("Using randomised SVD for whitening (faster, but less precise)")
  if(!fast_svd && .verbose)
    message("Using full SVD for whitening (slowed, but more precise)")
  c(X1, K) %<-% rs_prepare_whitening(
    x = target_mat,
    fast_svd = fast_svd,
    seed = random_seed,
    rank = nrow(target_mat),
    oversampling = NULL,
    n_power_iter = NULL
  )

  S7::prop(object, "processed_data")[["X1"]] <- X1
  S7::prop(object, "processed_data")[["K"]] <- K
  S7::prop(object, "params")["detection_method"] <- "ICA-based"

  return(object)
}

# component identification -----------------------------------------------------


ica_evaluate_comp <- S7::new_generic(
  name = "ica_evaluate_comp",
  dispatch_args = "object",
  fun = function() {

  }
)

S7::method(ica_evaluate_comp, bulk_coexp) <- function(object,
                                                      ica_type = c("logcosh", "exp"),
                                                      iter_params = list(
                                                        bootstrap = FALSE,
                                                        random_init = 50L,
                                                        folds = 10L,
                                                      ),
                                                      ncomp_params = list(
                                                        max_comp = 100L,
                                                        steps = 5L,
                                                        custom_seq = NULL
                                                      ),
                                                      ica_params = list(
                                                        maxit = 200L,
                                                        alpha = 1.0,
                                                        max_tol = 0.0001,
                                                        verbose = FALSE
                                                      ),
                                                      random_seed = 42L,
                                                      .verbose = TRUE) {
  checkmate::assertClass(object, "bixverse::bulk_coexp")
  checkmate::assertChoice(ica_fun, c("logcosh", "exp"))
  checkmate::qassert(.verbose, "B1")

  detection_method <- S7::prop(object, "params")[["detection_method"]]

  # Early return
  if (is.null(detection_method) &&
      detection_method != "ICA-based") {
    warning(
      paste(
        "This class does not seem to be set for ICA-based module detection",
        "Returning class as is."
      )
    )
    return(object)
  }

  # Prepare the n_comp vector
  n_comp_vector <- if (is.null(ncomp_params$custom_seq)) {
    with(ncomp_params, c(2, 3, 4, seq(
      from = 5, to = max_comp, by = steps
    )))
  } else {
    ncomp_params$custom_seq
  }
  if (.verbose)
    message(sprintf("Using a total of %i different n_comp parameters", length(n_comp_vector)))

  # Set up the loop
  if (iter_params$bootstrap) {
    no_ica_runs <- iter_params$random_init * iter_params$folds * length(n_comp_vector)
    if (.verbose)
      message(
        sprintf(
          "Using bootstrapping with %i folds and %i random initialisations for a total of %i ICA runs",
          iter_params$folds,
          iter_params$random_init,
          no_ica_runs
        )
      )
  } else {
    no_ica_runs <- iter_params$random_init * length(n_comp_vector)
    if (.verbose)
      message(
        sprintf(
          "Using %i random initialisations for a total of %i ICA runs",
          iter_params$random_init,
          no_ica_runs
        )
      )
  }

  all_scores <- c()
  all_convergence <- c()



}


# helpers ----------------------------------------------------------------------

## general ICA function --------------------------------------------------------

#' Fast ICA via Rust
#'
#' @description
#' This functions is a wrapper over the Rust implementation of fastICA. It has
#' the same two options `c("logcosh", "exp")` to run ICA in parallel modus. This
#' function expected
#'
#' @param X_norm Numeric matrix. Processed data. Output of
#' [bixverse::rs_prepare_whitening()].
#' @param K The K matrix. Whitening matrix. Output of
#' [bixverse::rs_prepare_whitening()].
#' @param n_icas Integer. Number of independent components to recover.
#' @param ica_fun String, element of `c("logcosh", "exp")`.
#' @param seed Integer. Seed to ensure reproducible results.
#' @param ica_params A list containing:
#' \itemize{
#'  \item maxit - Integer. Maximum number of iterations for ICA.
#'  \item alpha - Float. The alpha parameter for the logcosh version of ICA.
#'  Should be between 1 to 2.
#'  \item max_tol - Maximum tolerance of the algorithm.
#'  \item verbose - Controls verbosity of the function.
#' }
#' This list is optional. If a value cannot be found, default parameters will be
#' used.
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
                          ica_params = list(
                            maxit = 200L,
                            alpha = 1.0,
                            max_tol = 0.0001,
                            verbose = FALSE
                          )) {
  # Checks
  checkmate::assertMatrix(X_norm, mode = "numeric")
  checkmate::assertMatrix(K, mode = "numeric")
  checkmate::qassert(n_icas, sprintf("I1[2,%i]", ncol(X_norm)))
  checkmate::assertChoice(ica_fun, c("logcosh", "exp"))
  checkmate::qassert(seed, c("I1", "0"))
  # Function
  K <- matrix(K[1:n_icas, ], n_icas, dim(X_norm)[1])
  X1 <- K %*% X_norm
  set.seed(seed)
  w_init <- matrix(rnorm(n_icas * n_icas), nrow = n_icas, ncol = n_icas)

  c(a, converged) %<-% rs_fast_ica(
    X1,
    w_init,
    ica_type = ica_fun,
    ica_params = ica_params
  )

  w <- a %*% K
  S <- w %*% X_norm
  A <- t(w) %*% solve(w %*% t(w))

  res <- list(
    w = w,
    A = A,
    S = S,
    converged = converged
  )

  return(res)
}

