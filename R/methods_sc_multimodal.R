# multi modal functions --------------------------------------------------------

## adt -------------------------------------------------------------------------

### pca ------------------------------------------------------------------------

#' Calculate the PCA on top of the normalised ADT counts
#'
#' @description
#' This function will run PCA - via (randomised) SVD - on the normalised counts
#' and add the PCA results to the ScCache for the ADT counts.
#'
#' @param object `SingleCellsMultiModal` class with ADT counts added.
#' @param no_pcs Integer. Number of PCs to calculate.
#' @param features Optional string vector. If you want to subset to a specific
#' set of ADT probes (for example to exclude isotypes).
#' @param randomised_svd Boolean. Shall randomised SVD be used. Faster, but
#' less precise.
#' @param seed Integer. Controls reproducibility. Only relevant if
#' `randomised_svd = TRUE`.
#'
#' @returns The function will add the PCA factors, loadings and singular values
#' for the ADT data to the object.
#'
#' @export
calculate_pca_adt_sc <- S7::new_generic(
  name = "calculate_pca_adt_sc",
  dispatch_args = "object",
  fun = function(
    object,
    no_pcs,
    features = NULL,
    randomised_svd = FALSE,
    seed = 42L
  ) {
    S7::S7_dispatch()
  }
)

#' @method calculate_pca_adt_sc SingleCellsMultiModal
S7::method(calculate_pca_adt_sc, SingleCellsMultiModal) <- function(
  object,
  no_pcs,
  features = NULL,
  randomised_svd = FALSE,
  seed = 42L
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCellsMultiModal))
  checkmate::qassert(no_pcs, "I1")
  checkmate::qassert(features, c("0", "S+"))
  checkmate::qassert(randomised_svd, "B1")
  checkmate::qassert(seed, "I1")

  norm_counts <- get_sc_counts(
    object = object,
    assay = "norm",
    modality = "adt"
  )

  # force scaling
  zeallot::`%<-%`(
    c(pca_factors, pca_loadings, singular_values, is_scaled),
    if (randomised_svd) {
      rs_random_svd(
        norm_counts,
        scale = TRUE,
        rank = no_pcs,
        seed = seed,
        NULL,
        NULL
      )
    } else {
      rs_prcomp(x = norm_counts, scale = TRUE, top_pcs = no_pcs)
    }
  )

  object <- set_pca_factors(
    x = object,
    pca_factor = pca_factors,
    modality = "adt"
  )
  object <- set_pca_loadings(
    x = object,
    pca_loading = pca_loadings,
    modality = "adt"
  )
  object <- set_pca_singular_vals(
    x = object,
    singular_vals = singular_values[1:no_pcs],
    modality = "adt"
  )

  return(object)
}

## integration -----------------------------------------------------------------

### wnn ------------------------------------------------------------------------

#' Generate a weighted nearest neighbour (WNN) graph
#'
#' @description
#' This function implements the approach from Hao et al., to generate a weighted
#' nearest neighbour graph given two 'omics modalities (for example RNA and
#' ADT). Per-cell modality weights are computed by comparing within- and
#' cross-modality neighbourhood distances, and these weights are used to fuse
#' the two modalities into a single multimodal kNN graph.
#'
#' The per-modality kNN graphs are always recomputed internally from the chosen
#' embeddings at `knn_range` neighbours (see [params_sc_wnn()]); the kNN graphs
#' previously stored on the object via [bixverse::find_neighbours_sc()] are not
#' reused, as WNN requires a larger candidate pool than the default neighbour
#' search.
#'
#' The result is stored in the object's `other_data` under `"wnn"`, containing
#' the fused kNN graph as a `SingleCellNearestNeighbour` (with a kernelised
#' pseudo-distance metric), the sNN graph as an igraph and a table of per-cell
#' modality weights. This graph can subsequently be used for clustering and 2D
#' embedding via the `"wnn"` modality.
#'
#' @param object `SingleCellsMultiModal` class to which to add the WNN.
#' @param modality_1 String. First modality. Defaults to `"rna"`.
#' @param modality_2 String. Second modality. Defaults to `"adt"`.
#' @param embd_to_use_1 String. The embedding to use for the first modality.
#' Must be available in the object for `modality_1`. Defaults to `"pca"`.
#' @param embd_to_use_2 String. The embedding to use for the second modality.
#' Must be available in the object for `modality_2`. Defaults to `"pca"`.
#' @param no_embd_to_use_1 Optional integer. Number of embedding dimensions to
#' use for `embd_to_use_1`. If `NULL`, all will be used.
#' @param no_embd_to_use_2 Optional integer. Number of embedding dimensions to
#' use for `embd_to_use_2`. If `NULL`, all will be used.
#' @param wnn_params Named list. Controls the parameters for the WNN generation,
#' see [params_sc_wnn()].
#' @param full_snn Boolean. Shall the full shared nearest neighbour graph be
#' generated that generates edges between all cells instead of between only
#' neighbours.
#' @param pruning Numeric. Weights below this threshold will be set to 0 in the
#' generation of the sNN graph. Seurat uses for example 1/15 with k = 20.
#' @param snn_similarity String. One of `c("rank", "jaccard")`. The Jaccard
#' similarity calculates the Jaccard between the neighbours, whereas the rank
#' method calculates edge weights based on the ranking of shared neighbours. For
#' the rank method, the weight is determined by finding the shared neighbour
#' with the lowest combined rank across both cells, where lower-ranked (closer)
#' shared neighbours result in higher edge weights Both methods produce weights
#' normalised to the range `⁠[0, 1]`⁠.
#' @param seed Integer. For reproducibility.
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
#'
#' @returns The `SingleCellsMultiModal` object with the WNN graph and per-cell
#' modality weights added to `other_data[["wnn"]]`.
#'
#' @export
#'
#' @references Hao et al., Cell, 2021
generate_wnn_graph_sc <- S7::new_generic(
  name = "generate_wnn_graph_sc",
  dispatch_args = "object",
  fun = function(
    object,
    modality_1 = "rna",
    modality_2 = "adt",
    embd_to_use_1 = "pca",
    embd_to_use_2 = "pca",
    no_embd_to_use_1 = NULL,
    no_embd_to_use_2 = NULL,
    wnn_params = params_sc_wnn(),
    full_snn = TRUE,
    pruning = 1 / 15,
    snn_similarity = "jaccard",
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method generate_wnn_graph_sc SingleCellsMultiModal
S7::method(generate_wnn_graph_sc, SingleCellsMultiModal) <- function(
  object,
  modality_1 = "rna",
  modality_2 = "adt",
  embd_to_use_1 = "pca",
  embd_to_use_2 = "pca",
  no_embd_to_use_1 = NULL,
  no_embd_to_use_2 = NULL,
  wnn_params = params_sc_wnn(),
  full_snn = TRUE,
  pruning = 1 / 15,
  snn_similarity = "jaccard",
  seed = 42L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCellsMultiModal))
  checkmate::qassert(modality_1, "S1")
  checkmate::qassert(modality_2, "S1")
  checkmate::qassert(embd_to_use_1, "S1")
  checkmate::qassert(embd_to_use_2, "S1")
  checkmate::qassert(no_embd_to_use_1, c("I1", "0"))
  checkmate::qassert(no_embd_to_use_2, c("I1", "0"))
  assertScWnnParams(wnn_params)
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  # embeddings
  checkmate::assertTRUE(
    embd_to_use_1 %in% get_available_embeddings(object, modality = modality_1)
  )
  checkmate::assertTRUE(
    embd_to_use_2 %in% get_available_embeddings(object, modality = modality_2)
  )

  embd_1 <- get_embedding(
    x = object,
    embd_name = embd_to_use_1,
    modality = modality_1
  )
  embd_2 <- get_embedding(
    x = object,
    embd_name = embd_to_use_2,
    modality = modality_2
  )

  if (!is.null(no_embd_to_use_1)) {
    embd_1 <- embd_1[, 1:min(no_embd_to_use_1, ncol(embd_1))]
  }
  if (!is.null(no_embd_to_use_2)) {
    embd_2 <- embd_2[, 1:min(no_embd_to_use_2, ncol(embd_2))]
  }

  if (.verbose) {
    message("Running WNN graph generation.")
  }

  wnn_res <- rs_wnn(
    modality_emb_one = embd_1,
    modality_emb_two = embd_2,
    wnn_params = wnn_params,
    seed = seed,
    verbose = parse_verbosity(.verbose)
  )

  # store as a SingleCellNearestNeighbour
  used_cells <- get_cell_names(object, filtered = TRUE)
  wnn_knn <- new_sc_knn(
    knn_data = list(
      indices = wnn_res$indices,
      dist = wnn_res$dist,
      dist_metric = wnn_res$dist_metric
    ),
    used_cells = used_cells
  )

  wnn_weights <- data.table::data.table(
    cell_id = used_cells,
    weight_modality_1 = wnn_res$modality_one_weights,
    weight_modality_2 = wnn_res$modality_two_weights
  )

  if (.verbose) {
    message(sprintf(
      "Generating sNN graph (full: %s) from WNN graph.",
      full_snn
    ))
  }
  snn_graph_rs <- rs_sc_snn(
    knn_mat = get_knn_mat(wnn_knn),
    snn_method = snn_similarity,
    pruning = pruning,
    limited_graph = !full_snn,
    verbose = parse_verbosity(.verbose)
  )

  if (.verbose) {
    message("Transforming sNN data to igraph.")
  }
  snn_g <- igraph::make_empty_graph(n = nrow(embd_1), directed = FALSE)
  snn_g <- igraph::add_edges(
    snn_g,
    snn_graph_rs$edges,
    attr = list(weight = snn_graph_rs$weights)
  )

  other_data <- S7::prop(object, "other_data")
  other_data[["wnn"]] <- list(
    embeddings = list(),
    knn = wnn_knn,
    snn = snn_g,
    weights = wnn_weights,
    modality_1 = modality_1,
    modality_2 = modality_2
  )
  S7::prop(object, "other_data") <- other_data

  return(object)
}
