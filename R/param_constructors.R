# constructors -----------------------------------------------------------------

## miloR -----------------------------------------------------------------------

#' Wrapper function for parameters for MiloR
#'
#' @param prop Numeric. Proportion of cells to sample as neighbourhood indices.
#' Defaults to `0.2`. Must be in (0,1).
#' @param k_refine Integer. Number of neighbours to use for refinement.
#' Defaults to `20L`.
#' @param refinement_strategy String. Strategy for refining sampled indices.
#' One of `c("approximate", "bruteforce", "index")`. Defaults to
#' `"index"`.
#' @param index_type String. Type of kNN index to use. One of
#' `c("nndescent", "hnsw", "annoy", "exhaustive")`. Defaults to
#' `"nndescent"`. `"exhaustive"` scans every cell, so it returns the true
#' nearest neighbour rather than an approximation, at a cost that grows with
#' the number of cells.
#' @param knn List. Optional overrides for kNN parameters. See
#' [bixverse::params_knn_defaults()] for available parameters: `k`,
#' `knn_method`, `ann_dist`, `search_budget`, `n_trees`, `delta`,
#' `diversify_prob`, `ef_budget`, `m`, `ef_construction`, `ef_search`,
#' `n_list` and `n_probe`. `knn_method` is ignored in favour of `index_type`.
#'
#' @returns A list with the MiloR parameters and the full kNN block.
#'
#' @export
params_sc_miloR <- function(
  prop = 0.2,
  k_refine = 20L,
  refinement_strategy = c("index", "approximate", "bruteforce"),
  index_type = c("nndescent", "hnsw", "annoy", "exhaustive"),
  knn = list()
) {
  refinement_strategy <- match.arg(refinement_strategy)
  index_type <- match.arg(index_type)
  checkmate::qassert(prop, "N1(0,1)")
  checkmate::qassert(k_refine, "I1[1,)")
  checkmate::assertList(knn, names = "unique")
  checkmate::assertSubset(names(knn), names(params_knn_defaults()))

  knn_params <- modifyList(
    params_knn_defaults(),
    knn,
    keep.null = TRUE
  )
  assertKnnParams(knn_params)

  c(
    list(
      prop = prop,
      k_refine = k_refine,
      refinement_strategy = refinement_strategy,
      index_type = index_type
    ),
    knn_params
  )
}
