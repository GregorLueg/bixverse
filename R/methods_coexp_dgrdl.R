# grid search ------------------------------------------------------------------

#' @title Grid search over DGRDL parameters
#'
#' @description
#' This function allows you to quickly iterate over different initial seeds,
#' number of neighbours for the KNN graph and dictionary sizes to identify
#' optimal hyperparameters for your DGRDL problem.
#'
#' @param object The class, see [bixverse::bulk_coexp()]. Ideally, you
#' should run [bixverse::preprocess_bulk_coexp()] before applying this function.
#' @param neighbours_vec Integer vector. The different k nearest neighbours
#' to test.
#' @param dict_size_vec Integer vector. The different dictionary sizes to test
#' for.
#' @param seed_vec Integer vector. The different initial seeds to test for the
#' dictionary generation.
#' @param dgrdl_params List. Output of [bixverse::params_dgrdl()]:
#' \itemise {
#'   \item sparsity - Integer. Sparsity constraint (max non-zero coefficients
#'   per signal)
#'   \item dict size - Integer. Will be ignored by this function and the
#'   `dict_size_vec` vector will be used.
#'   \item alpha - Float. Sample context regularisation weight.
#'   \item beta - Float. Feature effect regularisation weight.
#'   \item max_iter - Integer. Maximum number of iterations for the main
#'   algorithm.
#'   \item k_neighbours - Integer. Will be ignored by this function and the
#'   `neighbours_vec` will be used.
#'   \item admm_iter - Integer. ADMM iterations for sparse coding.
#'   \rho rho - Float. ADMM step size.
#' }
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @return `bulk_coexp` with the grid search results added to the class.
#'
#' @export
dgrdl_grid_search <- S7::new_generic(
  name = "dgrdl_grid_search",
  dispatch_args = "object",
  fun = function(
    object,
    neighbours_vec,
    dict_size_vec,
    seed_vec,
    dgrdl_params = params_dgrdl(),
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @export
#'
#' @import data.table
#' @importFrom magrittr `%>%`
#'
#' @method dgrdl_grid_search bulk_coexp
S7::method(dgrdl_grid_search, bulk_coexp) <- function(
  object,
  neighbours_vec,
  dict_size_vec,
  seed_vec,
  dgrdl_params = params_dgrdl(),
  .verbose = TRUE
) {
  # checks
  checkmate::assertClass(object, "bixverse::bulk_coexp")
  checkmate::qassert(neighbours_vec, "I+")
  checkmate::qassert(dict_size_vec, "I+")
  checkmate::qassert(seed_vec, "I+")
  assertDGRDLparams(dgrdl_params)
  checkmate::qassert(.verbose, "B1")

  # function body
  if (purrr::is_empty(S7::prop(object, "processed_data")[["processed_data"]])) {
    warning("No pre-processed data found. Defaulting to the raw data.")
    target_mat <- S7::prop(object, "raw_data")
  } else {
    target_mat <- S7::prop(object, "processed_data")[["processed_data"]]
  }

  total_params <- length(neighbours_vec) *
    length(dict_size_vec) *
    length(seed_vec)

  if (.verbose) {
    message(
      "A total of %i parameters will be tested in a grid search for DGRDL.",
      total_params
    )
  }

  # do the grid search
  grid_search_res <- rs_sparse_dict_dgrdl_grid_search(
    x = target_mat,
    dgrdl_params = dgrdl_params,
    seeds = seed_vec,
    dict_sizes = dict_size_vec,
    k_neighbours_vec = neighbours_vec,
    verbose = .verbose
  ) %>%
    data.table::setDT()

  grid_search_params <- with(
    dgrdl_params,
    list(
      sparsity = sparsity,
      tested_dict_sizes = dict_size_vec,
      alpha = alpha,
      beta = beta,
      max_iters = max_iters,
      tested_k_neighbours = neighbours_vec,
      admm_iter = admm_iter,
      rho = rho,
      tested_seeds = seed_vec
    )
  )

  S7::prop(object, "params")[["grid_search_params"]] <- grid_search_params
  S7::prop(object, "outputs")[["grid_search_res"]] <- grid_search_res

  return(object)
}

## plotting --------------------------------------------------------------------

# dgrdl ------------------------------------------------------------------------

# specific getters -------------------------------------------------------------

#' @title Get the grid search results
#'
#' @description
#' Getter function to extract the grid search results. If not found will return
#' `NULL`.
#'
#' @param object The class, see [bixverse::bulk_coexp()].
#'
#' @return data.table with the grid search results (if found. Otherwise `NULL`.)
#'
#' @export
get_grid_search_res <- S7::new_generic(
  name = "get_grid_search_res",
  dispatch_args = "object",
  fun = function(object) {
    S7::S7_dispatch()
  }
)

#' @export
#'
#' @method get_grid_search_res bulk_coexp
S7::method(get_grid_search_res, bulk_coexp) <- function(object) {
  checkmate::assertClass(object, "bixverse::bulk_coexp")

  grid_search_res <- S7::prop(object, "outputs")[["grid_search_res"]]

  if (is.null(grid_search_res)) {
    warning(paste(
      "No grid search results found.",
      "Did you run dgrdl_grid_search()? Returning NULL."
    ))
  }

  return(grid_search_res)
}
