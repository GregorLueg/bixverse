# single cell aggregation methods ----------------------------------------------

## meta cells ------------------------------------------------------------------

#' Generate meta cells and return the matrices
#'
#' @description
#' This function implements the meta cell aggregation from Morabito, et al.
#' The generation of metacells is a useful approach for subsequent application
#' of for example correlation-based methods to identify co-regulated genes.
#'
#' @param object `single_cell_exp` class.
#' @param sc_meta_cell_params List. Output of [bixverse::params_sc_metacells()].
#' A list with the following items:
#' \itemize{
#'   \item max_shared - Maximum number of allowed shared neighbours for
#'   the meta cell to be considered.
#'   \item target_no_metacells - Number of target meta cells you wish to reach.
#'   \item max_iter - Maximum number of iterations you want to use for the
#'   algorithm.
#'   \item k - Number of neighbours for the kNN search. Only relevant if you
#'   set regenerate_knn to `TRUE`.
#'   \item n_trees - Integer. Number of trees to use for the annoy algorithm.
#'   Only relevant if you set regenerate_knn to `TRUE`.
#'   \item search_budget - Integer. Search budget per tree for the annoy
#'   algorithm. Only relevant if you set regenerate_knn to `TRUE`.
#'   \item knn_algorithm - String. Which of the two implemented kNN algorithms
#'   to use. Defaults to `"annoy"`. Only relevant if you set regenerate_knn to
#'   `TRUE`.
#' }
#' @param embd_to_use String. The embedding to use. Atm, the only option is
#' `"pca"`. Only relevant if you set regenerate_knn to `TRUE`.
#' @param no_embd_to_use Optional integer. Number of embedding dimensions to
#' use. If `NULL` all will be used. Only relevant if you set regenerate_knn to
#' `TRUE`.
#' @param target_size Numeric. The library target size to normalise the meta
#' cells to.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @references
#' Morabito, et al. Cell Rep Methods, 2023
get_meta_cells_sc <- S7::new_generic(
  name = "get_meta_cells_sc",
  dispatch_args = "object",
  fun = function(
    object,
    sc_meta_cell_params = params_sc_metacells(),
    regenerate_knn = FALSE,
    embd_to_use = "pca",
    no_embd_to_use = NULL,
    target_size = 1e4,
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method get_meta_cells_sc single_cell_exp
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(get_meta_cells_sc, single_cell_exp) <- function(
  object,
  sc_meta_cell_params = params_sc_metacells(),
  regenerate_knn = FALSE,
  embd_to_use = "pca",
  no_embd_to_use = NULL,
  target_size = 1e4,
  seed = 42L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertClass(object, "bixverse::single_cell_exp")
  assertScMetacells(sc_meta_cell_params)
  checkmate::qassert(regenerate_knn, "B1")
  checkmate::assertChoice(embd_to_use, c("pca"))
  checkmate::qassert(no_embd_to_use, c("I1", "0"))
  checkmate::qassert(target_size, "N1")
  checkmate::qassert(.verbose, "B1")

  # if the kNN graph shall be regenerated, get the emedding here...
  if (regenerate_knn) {
    embd <- switch(embd_to_use, pca = get_pca_factors(object))
    # early return
    if (is.null(embd)) {
      warning(
        paste(
          "The desired embedding was not found. Please check the parameters.",
          "Returning NULL."
        )
      )

      return(NULL)
    }

    if (!is.null(no_embd_to_use)) {
      to_take <- min(c(no_embd_to_use, ncol(embd)))
      embd <- embd[, 1:to_take]
    }
    knn_data <- NULL
  } else {
    embd <- NULL
    knn_data <- get_knn_mat(object)

    if (is.null(knn_data)) {
      warning(
        paste(
          "No kNN data could be found on the object. Set regenerate_knn to",
          "TRUE or generate the kNN matrix via other means",
          "Returning NULL."
        )
      )
      return(NULL)
    }
  }

  meta_cell_data <- rs_get_metacells(
    f_path = get_rust_count_cell_f_path(object),
    knn_mat = knn_data,
    embd = embd,
    meta_cell_params = sc_meta_cell_params,
    target_size = target_size,
    seed = seed,
    verbose = .verbose
  )

  c(raw_counts, norm_counts) %<-% get_meta_cell_matrices(meta_cell_data)

  no_gen_meta_cells <- nrow(raw_counts)

  colnames(raw_counts) <- colnames(norm_counts) <- get_gene_names(object)
  n_digits <- nchar(as.character(no_gen_meta_cells))
  format_str <- sprintf("meta_cell_%%0%dd", n_digits)
  rownames(raw_counts) <- rownames(norm_counts) <- sprintf(
    format_str,
    1:no_gen_meta_cells
  )

  res <- list(meta_raw_counts = raw_counts, meta_norm_counts = norm_counts)

  res
}
