# methods related to batch corrections -----------------------------------------

## kBET scores -----------------------------------------------------------------

#' Calculate kBET scores
#'
#' @description
#' This function calculates the k-nearest neighbour batch-effect test (kBET).
#' Briefly, the function leverages a Chi-Square statistic to calculate the
#' differences in batch proportions observed in the neighbourhood of a given
#' cell with the overall batch proportions. If the test is significant for that
#' cell it indicates poor mixing for that cell specifically.
#' Large number of positive tests indicate bad mixing overall. For more details,
#' please see Büttner et al.
#'
#' @param object `single_cell_exp` class.
#' @param batch_column String. The column with the batch information in the
#' obs data of the class.
#' @param threshold Numeric. Number between 0 and 1. Below this threshold, the
#' test is considered significant. Defaults to `0.05`.
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @returns A list with the following elements
#' \itemize{
#'   \item kbet_score - Number of significant tests over all cells. 0 indicates
#'   perfect mixing, 1 indicates basically zero mixing between batches.
#'   \item significant_tests - Boolean indicating for which cells the statistic
#'   was below the threshold
#'   \item chisquare_pvals - The p-values of the ChiSquare test.
#' }
#'
#' @export
#'
#' @references Büttner, et al., Nat. Methods, 2019
calculate_kbet_sc <- S7::new_generic(
  name = "calculate_kbet_sc",
  dispatch_args = "object",
  fun = function(
    object,
    batch_column,
    threshold = 0.05,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)


#' @method calculate_kbet_sc single_cell_exp
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(calculate_kbet_sc, single_cell_exp) <- function(
  object,
  batch_column,
  threshold = 0.05,
  .verbose = TRUE
) {
  # check
  checkmate::assertTRUE(S7::S7_inherits(object, single_cell_exp))
  checkmate::qassert(batch_column, "S1")
  checkmate::qassert(threshold, "N1[0, 1]")
  checkmate::qassert(.verbose, "B1")

  # function
  batch_index <- unlist(object[[batch_column]])

  if (!length(levels(factor(batch_index))) > 1) {
    warning("The batch column only has one batch. Returning NULL")
    return(NULL)
  }

  knn_mat <- get_knn_mat(object)

  if (is.null(knn_mat)) {
    warning("No kNN matrix found to calculate kBET. Returning NULL")
    return(NULL)
  }

  rs_res <- rs_kbet(
    knn_mat = knn_mat,
    batch_vector = as.integer(factor(batch_index))
  )

  res <- list(
    kbet_score = sum(rs_res <= threshold) / length(rs_res),
    significant_tests = rs_res <= threshold,
    chisquare_pvals = rs_res
  )

  return(res)
}

## BBKNN -----------------------------------------------------------------------

#' Run BBKNN
#'
#' @description
#' This function implements the batch-balanced k-nearest neighbour algorithm
#' from Polański, et al. Briefly, the algorithm generate a KNN index on a per
#' batch basis and identifies the neighbours of cells for each individual index.
#' Subsequently, it leverages UMAP connectivity calculations to reduce spurious
#' connections. For more details, please refer to Polański, et al.
#'
#' @param object `single_cell_exp` class.
#' @param batch_column String. The column with the batch information in the
#' obs data of the class.
#' @param no_neighbours_to_keep Integer. Maximum number of neighbours to keep
#' from the BBKNN algorithm. Due to generating neighbours for each batch, there
#' might be a large number of generated neighbours. This will only keep the top
#' `no_neighbours_to_keep` neighbours.
#' @param embd_to_use String. The embedding to use. Atm, the only option is
#' `"pca"`.
#' @param no_embd_to_use Optional integer. Number of embedding dimensions to
#' use. If `NULL` all will be used.
#' @param sc_bbknn_params A list, please see [bixverse::params_sc_bbknn()]. The
#' list has the following parameters:
#' \itemize{
#'   \item neighbours_within_batch - Integer. Number of neighbours to consider
#'   per batch.
#'   \item knn_method - String. One of `c("annoy", "hnsw")`. Defaults to
#'   `"annoy"`.
#'   \item ann_dist - String. One of `c("cosine", "euclidean")`. The distance
#'   metric to be used for the approximate neighbour search. Defaults to
#'   `"cosine"`.
#'   \item set_op_mix_ratio - Numeric. Mixing ratio between union (1.0) and
#'   intersection (0.0).
#'   \item local_connectivity - Numeric. UMAP connectivity computation
#'   parameter, how many nearest neighbours of each cell are assumed to be fully
#'   connected.
#'   \item annoy_n_trees - Integer. Number of trees to use in the generation of
#'   the Annoy index.
#'   \item search_budget - Integer. Search budget per tree for the `annoy`
#'   algorithm.
#'   \item trim - Optional integer. Trim the neighbours of each cell to these
#'   many top connectivities. May help with population independence and improve
#'   the tidiness of clustering. If `NULL`, it defaults to
#'   `10 * neighbours_within_batch`.
#' }
#' @param seed Integer. Random seed.
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @returns The object with added kNN matrix based on BBKNN and the graph based
#' on the returned connectivities of the algorithm.
#'
#' @export
#'
#' @references Polański, et al., Bioinformatics, 2020
bbknn_sc <- S7::new_generic(
  name = "bbknn_sc",
  dispatch_args = "object",
  fun = function(
    object,
    batch_column,
    no_neighbours_to_keep = 15L,
    embd_to_use = "pca",
    no_embd_to_use = NULL,
    sc_bbknn_params = params_sc_bbknn(),
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method bbknn_sc single_cell_exp
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(bbknn_sc, single_cell_exp) <- function(
  object,
  batch_column,
  no_neighbours_to_keep = 15L,
  embd_to_use = "pca",
  no_embd_to_use = NULL,
  sc_bbknn_params = params_sc_bbknn(),
  seed = 42L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, single_cell_exp))
  checkmate::qassert(batch_column, "S1")
  checkmate::assertChoice(embd_to_use, c("pca"))
  checkmate::qassert(no_embd_to_use, c("I1", "0"))
  assertScBbknn(sc_bbknn_params)
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, "B1")

  # function body
  if (!is.null(get_knn_mat(object))) {
    warning("Prior kNN matrix found. Will be overwritten.")
  }

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

  batch_index <- unlist(object[[batch_column]])

  if (!length(levels(factor(batch_index))) > 1) {
    warning("The batch column only has one batch. Returning object as is.")
    return(object)
  }

  no_generated_neighbours <- length(levels(factor(batch_index))) *
    sc_bbknn_params$neighbours_within_batch

  if (no_neighbours_to_keep > no_generated_neighbours) {
    warning(paste(
      "The number of desired neighbours cannot be generated with these BBKNN",
      "parameters (too few generated neighbours).",
      "Please adopt neighbours_within_batch accordingly.",
      "Returning all neighbours from BBKNN."
    ))
  }

  if (.verbose) {
    message("Running BBKNN algorithm.")
  }

  bbknn_res <- rs_bbknn(
    embd = embd,
    batch_labels = as.integer(batch_index),
    bbknn_params = sc_bbknn_params,
    seed = seed,
    verbose = .verbose
  )

  knn_mat <- if (no_generated_neighbours <= no_neighbours_to_keep) {
    matrix(data = bbknn_res$distances$indices, nrow = nrow(embd), byrow = TRUE)
  } else {
    rs_bbknn_filtering(
      indptr = bbknn_res$distances$indptr,
      indices = bbknn_res$distances$indices,
      no_neighbours_to_keep = no_neighbours_to_keep
    )
  }

  storage.mode(knn_mat) <- "integer"

  object <- set_knn(object, knn_mat = knn_mat)

  if (.verbose) {
    message(paste(
      "Generating graph based on BBKNN connectivities.",
      "Weights will be based on the connectivities and not shared nearest neighour",
      "calculations."
    ))
  }

  sparse_mat <- Matrix::sparseMatrix(
    i = rep(
      seq_along(bbknn_res$connectivities$indptr[-1]),
      diff(bbknn_res$connectivities$indptr)
    ),
    j = bbknn_res$connectivities$indices + 1,
    x = bbknn_res$connectivities$data,
    dims = c(bbknn_res$connectivities$nrow, bbknn_res$connectivities$ncol),
    index1 = TRUE
  )

  snn_graph <- igraph::graph_from_adjacency_matrix(
    sparse_mat,
    mode = "max",
    weighted = TRUE
  )

  object <- set_snn_graph(object, snn_graph = snn_graph)

  return(object)
}
