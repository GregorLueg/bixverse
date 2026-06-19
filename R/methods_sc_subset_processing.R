## hvg -------------------------------------------------------------------------

# generic found in R/base_generics_sc.R

#' @method find_hvg_sc SingleCellsSubset
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(find_hvg_sc, SingleCellsSubset) <- function(
  object,
  hvg_no = 2000L,
  hvg_params = params_sc_hvg(),
  streaming = NULL,
  .verbose = TRUE
) {
  checkmate::assertClass(object, "bixverse::SingleCellsSubset")
  checkmate::qassert(hvg_no, "I1")
  assertScHvg(hvg_params)
  checkmate::qassert(streaming, c("B1", "0"))
  checkmate::qassert(.verbose, c("B1", "I1[0, 2]"))

  streaming <- auto_streaming(
    n_cells = nrow(object),
    streaming = streaming,
    .verbose = .verbose
  )

  res <- with(
    hvg_params,
    rs_sc_hvg(
      f_path_gene = get_rust_count_gene_f_path(object),
      hvg_method = method,
      cell_indices = get_cells_to_keep(object),
      loess_span = loess_span,
      n_bins = num_bin,
      binning = bin_method,
      clip_max = NULL,
      streaming = streaming,
      verbose = parse_verbosity(.verbose)
    )
  )

  hvg <- switch(
    hvg_params$method,
    "vst" = order(res$var_std, decreasing = TRUE)[1:hvg_no],
    "dispersion" = order(res$dispersion, decreasing = TRUE)[1:hvg_no],
    "meanvarbin" = order(res$dispersion_scaled, decreasing = TRUE)[1:hvg_no],
    stop("Unknown HVG method: ", hvg_params$method)
  )

  object <- set_hvg(object, hvg = hvg)

  return(object)
}

## dimension reduction and knn/snn ---------------------------------------------

### pca ------------------------------------------------------------------------

# generic found in R/base_generics_sc.R

#' @method calculate_pca_sc SingleCellsSubset
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(calculate_pca_sc, SingleCellsSubset) <- function(
  object,
  no_pcs,
  randomised_svd = TRUE,
  sparse_svd = FALSE,
  hvg = NULL,
  seed = 42L,
  .verbose = TRUE
) {
  checkmate::assertClass(object, "bixverse::SingleCellsSubset")
  checkmate::qassert(no_pcs, "I1")
  checkmate::qassert(randomised_svd, "B1")
  checkmate::qassert(sparse_svd, "B1")
  checkmate::qassert(hvg, c("I+", "0"))
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  if ((length(get_hvg(object)) == 0) && is.null(hvg)) {
    warning(paste(
      "No HVGs identified in the object nor provided.",
      "Please run find_hvg_sc() or provide the indices of the HVG",
      "Returning object as is."
    ))
    return(object)
  }

  selected_hvg <- if (!is.null(hvg)) {
    if (.verbose) {
      message(
        paste(
          "HVGs provided.",
          "Will use these ones and set the internal HVG to the provided genes."
        )
      )
    }
    object <- set_hvg(object, hvg) # this one deals with zero/one indexing internally
    hvg - 1L
  } else {
    get_hvg(object)
  }

  # swap to sparse SVD for large data sets
  n_cells <- length(get_cells_to_keep(object))

  if (n_cells > 500000 & !sparse_svd) {
    message(paste(
      "More than 500,000 cells with sparse SVD = FALSE",
      "Setting sparse SVD to TRUE to avoid high memory pressure."
    ))

    sparse_svd <- TRUE
  }

  # dense path
  if (!sparse_svd) {
    if (.verbose) {
      message(
        sprintf(
          "Using dense SVD solving on scaled data on %i HVG.",
          length(selected_hvg)
        )
      )
    }
    zeallot::`%<-%`(
      c(pca_factors, pca_loadings, singular_values, scaled),
      rs_sc_pca(
        f_path_gene = get_rust_count_gene_f_path(object),
        no_pcs = no_pcs,
        random_svd = randomised_svd,
        cell_indices = get_cells_to_keep(object),
        gene_indices = selected_hvg,
        seed = seed,
        return_scaled = FALSE,
        verbose = parse_verbosity(.verbose)
      )
    )

    object <- set_pca_factors(object, pca_factors)
    object <- set_pca_loadings(object, pca_loadings)
    object <- set_pca_singular_vals(object, singular_values[1:no_pcs])

    return(object)
  } else {
    if (.verbose) {
      message(
        sprintf(
          "Using sparse SVD solving on scaled data on %i HVG.",
          length(selected_hvg)
        )
      )
    }

    zeallot::`%<-%`(
      c(sparse_pca_factors, sparse_pca_loadings, sparse_pca_eigenvals),
      rs_sc_pca_sparse(
        f_path_gene = bixverse:::get_rust_count_gene_f_path(object),
        no_pcs = no_pcs,
        random_svd = randomised_svd,
        cell_indices = get_cells_to_keep(object),
        gene_indices = selected_hvg,
        seed = seed,
        verbose = parse_verbosity(.verbose)
      )
    )

    object <- set_pca_factors(object, sparse_pca_factors)
    object <- set_pca_loadings(object, sparse_pca_loadings)
    object <- set_pca_singular_vals(
      object,
      sparse_pca_eigenvals[1:no_pcs]
    )

    return(object)
  }
}

### neighbours -----------------------------------------------------------------

# generic found in R/base_generics_sc.R
# method shared across SingleCellsSubset and MetaCells

S7::method(find_neighbours_sc, ScOrMc) <- function(
  object,
  embd_to_use = "pca",
  no_embd_to_use = NULL,
  modality = c("rna", "adt"),
  neighbours_params = params_sc_neighbours(),
  seed = 42L,
  .verbose = TRUE
) {
  modality <- match.arg(modality)

  checkmate::assertTRUE(
    S7::S7_inherits(object, SingleCellsSubset) ||
      S7::S7_inherits(object, MetaCells)
  )
  checkmate::qassert(embd_to_use, "S1")
  checkmate::qassert(no_embd_to_use, c("I1", "0"))
  assertScNeighbours(neighbours_params)
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  if (
    modality != "rna" && !S7::S7_inherits(object, SingleCellsSubsetMultiModal)
  ) {
    stop(sprintf(
      "modality = '%s' is only supported for SingleCellsSubsetMultiModal.",
      modality
    ))
  }

  if (!embd_to_use %in% get_available_embeddings(object, modality = modality)) {
    warning("The desired embedding was not found. Returning class as is.")
    return(object)
  }

  embd <- get_embedding(
    x = object,
    embd_name = embd_to_use,
    modality = modality
  )
  if (!is.null(no_embd_to_use)) {
    embd <- embd[, 1:min(no_embd_to_use, ncol(embd))]
  }

  if (.verbose) {
    message(sprintf(
      "Generating kNN data with %s method.",
      neighbours_params$knn_algorithm
    ))
  }
  knn_data <- generate_sc_knn(
    data = embd,
    neighbours_params = neighbours_params,
    seed = seed,
    .verbose = .verbose
  )
  object <- set_knn(object, knn_data, modality = modality)

  if (.verbose) {
    message(sprintf(
      "Generating sNN graph (full: %s).",
      neighbours_params$full_snn
    ))
  }
  snn_graph_rs <- with(
    neighbours_params,
    rs_sc_snn(
      knn_mat = get_knn_mat(knn_data),
      snn_method = snn_similarity,
      pruning = pruning,
      limited_graph = !full_snn,
      verbose = parse_verbosity(.verbose)
    )
  )

  if (.verbose) {
    message("Transforming sNN data to igraph.")
  }
  snn_g <- igraph::make_empty_graph(n = nrow(embd), directed = FALSE)
  snn_g <- igraph::add_edges(
    snn_g,
    snn_graph_rs$edges,
    attr = list(weight = snn_graph_rs$weights)
  )

  object <- set_snn_graph(object, snn_graph = snn_g, modality = modality)

  return(object)
}

### clustering -----------------------------------------------------------------

# generic found in R/base_generics_sc.R
# method shared across SingleCellsSubset and MetaCells

S7::method(find_clusters_sc, ScOrMc) <- function(
  object,
  cluster_algorithm = c("leiden", "louvain"),
  res = 1.0,
  name = "leiden_clustering",
  modality = c("rna", "adt", "wnn"),
  seed = 42L
) {
  cluster_algorithm <- match.arg(cluster_algorithm)
  modality <- match.arg(modality)

  checkmate::assertTRUE(
    S7::S7_inherits(object, SingleCellsSubset) ||
      S7::S7_inherits(object, MetaCells)
  )
  checkmate::qassert(res, "N1")
  checkmate::qassert(name, "S1")
  checkmate::assertChoice(cluster_algorithm, c("leiden", "louvain"))
  checkmate::qassert(seed, "I1")

  if (
    modality != "rna" && !S7::S7_inherits(object, SingleCellsSubsetMultiModal)
  ) {
    stop(sprintf(
      "modality = '%s' is only supported for SingleCellsSubsetMultiModal.",
      modality
    ))
  }

  snn_graph <- get_snn_graph(object, modality = modality)
  if (is.null(snn_graph)) {
    warning(
      paste(
        "No sNN graph found. Did you run find_neighbours_sc().",
        "Returning class as is."
      )
    )
    return(object)
  }

  # set a global seed to ensure more reproducibility here
  set.seed(seed)
  clusters <- switch(
    cluster_algorithm,
    leiden = igraph::cluster_leiden(
      graph = snn_graph,
      objective_function = "modularity",
      resolution = res
    ),
    louvain = igraph::cluster_louvain(graph = snn_graph, resolution = res)
  )

  object[[name]] <- .remap_communities_by_size(clusters$membership)

  object
}

### fast clustering ------------------------------------------------------------

#' Run fast Louvain clustering on a SingleCellsSubset object
#'
#' @description
#' Runs k-means on the chosen embedding, builds a kNN graph on the centroids,
#' applies Louvain clustering and propagates memberships back to the cells.
#' Optionally runs a grid over multiple seeds and returns stability statistics.
#'
#' @param object `SingleCellsSubset` class.
#' @param embd_to_use String. Embedding name. Defaults to `"pca"`.
#' @param no_embd_to_use Optional integer. Number of dimensions to keep.
#' @param resolutions Numeric vector. Louvain resolutions.
#' @param km_type String. One of `c("kmeans", "minibatch")`. The former runs
#' standard k-means, the latter a mini-batch version that can be useful for
#' large data sets.
#' @param n_centroids Optional integer. Number of k-means centroids. Defaults
#' to `sqrt(n_cells)` Rust-side if `NULL`.
#' @param fc_params List. Output of [params_sc_fast_cluster()].
#' @param snn Boolean. Convert kNN to sNN.
#' @param return_kmeans Boolean. Return k-means assignments and centroids.
#' @param grid_search Boolean. Run multi-seed grid version.
#' @param no_seeds Integer. Number of additional seeds (only used when
#' `grid_search = TRUE`).
#' @param seed Integer. Reproducibility.
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
#'
#' @returns `SingleCellFastClusters` S3 object with:
#' \describe{
#'   \item{memberships}{data.table with `cell_idx` and one column per
#'   resolution (`res_<value>`).}
#'   \item{stats}{data.table of grid statistics, or `NULL`.}
#'   \item{k_means_cluster}{Integer vector of k-means assignments, or `NULL`.}
#'   \item{centroids}{Numeric matrix of centroids, or `NULL`.}
#'   \item{resolutions}{Resolutions used.}
#' }
#' with `cell_indices` stored as an attribute (0-indexed).
#'
#' @export
fast_cluster_sc <- S7::new_generic(
  name = "fast_cluster_sc",
  dispatch_args = "object",
  fun = function(
    object,
    embd_to_use = "pca",
    no_embd_to_use = NULL,
    resolutions = c(2.0, 1.0, 0.5),
    km_type = c("kmeans", "minibatch"),
    n_centroids = NULL,
    fc_params = params_sc_fast_cluster(),
    snn = TRUE,
    return_kmeans = FALSE,
    grid_search = FALSE,
    no_seeds = 10L,
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method fast_cluster_sc SingleCellsSubset
#'
#' @export
S7::method(fast_cluster_sc, SingleCellsSubset) <- function(
  object,
  embd_to_use = "pca",
  no_embd_to_use = NULL,
  resolutions = c(2.0, 1.0, 0.5),
  km_type = c("kmeans", "minibatch"),
  n_centroids = NULL,
  fc_params = params_sc_fast_cluster(),
  snn = TRUE,
  return_kmeans = FALSE,
  grid_search = FALSE,
  no_seeds = 10L,
  seed = 42L,
  .verbose = TRUE
) {
  km_type <- match.arg(km_type)

  checkmate::assertTRUE(S7::S7_inherits(object, SingleCellsSubset))
  checkmate::qassert(embd_to_use, "S1")
  checkmate::qassert(no_embd_to_use, c("I1", "0"))
  checkmate::qassert(resolutions, "N+")
  checkmate::qassert(n_centroids, c("I1", "0"))
  checkmate::qassert(snn, "B1")
  checkmate::qassert(return_kmeans, "B1")
  checkmate::qassert(grid_search, "B1")
  checkmate::qassert(no_seeds, "I1")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  if (!embd_to_use %in% get_available_embeddings(object)) {
    stop(sprintf("Embedding '%s' was not found.", embd_to_use))
  }

  embd <- get_embedding(x = object, embd_name = embd_to_use)
  if (!is.null(no_embd_to_use)) {
    embd <- embd[, 1:min(no_embd_to_use, ncol(embd))]
  }

  cells_to_use <- get_cells_to_keep(object)

  if (grid_search) {
    res <- rs_fast_cluster_sc_grid(
      embd = embd,
      km_type = km_type,
      resolutions = resolutions,
      n_centroids = n_centroids,
      fc_params = fc_params,
      snn = snn,
      return_kmeans = return_kmeans,
      no_seeds = no_seeds,
      seed = seed,
      verbose = parse_verbosity(.verbose)
    )
    memberships <- res$membership$memberships
    stats <- data.table::as.data.table(res$membership$stats)
    stats[, resolution := resolutions]
    data.table::setcolorder(stats, "resolution")
  } else {
    res <- rs_fast_cluster_sc(
      embd = embd,
      km_type = km_type,
      resolutions = resolutions,
      n_centroids = n_centroids,
      fc_params = fc_params,
      snn = snn,
      return_kmeans = return_kmeans,
      seed = seed,
      verbose = parse_verbosity(.verbose)
    )
    memberships <- res$membership
    stats <- NULL
  }

  membership_dt <- data.table::data.table(cell_idx = cells_to_use + 1L)
  membership_dt[, (paste0("res_", resolutions)) := memberships]

  structure(
    list(
      memberships = membership_dt,
      stats = stats,
      k_means_cluster = res$k_means_cluster,
      centroids = res$centroids,
      resolutions = resolutions
    ),
    cell_indices = cells_to_use,
    class = "SingleCellFastClusters"
  )
}
