# single cell subset processing methods ----------------------------------------

## hvg -------------------------------------------------------------------------

# generic in R/base_generics_sc.R

#' @method find_hvg_sc SingleCellsSubset
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
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

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

  set_hvg(object, hvg = hvg)
}

## pca -------------------------------------------------------------------------

# generic in R/base_generics_sc.R

#' @method calculate_pca_sc SingleCellsSubset
#'
#' @importFrom zeallot `%<-%`
S7::method(calculate_pca_sc, SingleCellsSubset) <- function(
  object,
  no_pcs,
  pca_params = params_sc_pca(),
  sparse_svd = FALSE,
  hvg = NULL,
  seed = 42L,
  .verbose = TRUE
) {
  checkmate::assertClass(object, "bixverse::SingleCellsSubset")
  checkmate::qassert(no_pcs, "I1")
  assertScPca(pca_params)
  checkmate::qassert(sparse_svd, "B1")
  checkmate::qassert(hvg, c("I+", "0"))
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  if ((length(get_hvg(object)) == 0) && is.null(hvg)) {
    warning(paste(
      "No HVGs identified in the object nor provided.",
      "Please run find_hvg_sc() or provide the indices of the HVG.",
      "Returning object as is."
    ))
    return(object)
  }

  selected_hvg <- if (!is.null(hvg)) {
    if (.verbose) {
      message(paste(
        "HVGs provided.",
        "Will use these and overwrite the internal HVG."
      ))
    }
    object <- set_hvg(object, hvg) # set_hvg takes 1-indexed, stores 0-indexed
    hvg - 1L
  } else {
    get_hvg(object) # already 0-indexed
  }

  n_cells <- length(get_cells_to_keep(object))
  if (n_cells > 500000L && !sparse_svd) {
    message(paste(
      "More than 500,000 cells with sparse SVD = FALSE.",
      "Setting sparse SVD to TRUE to avoid high memory pressure."
    ))
    sparse_svd <- TRUE
  }

  if (!sparse_svd) {
    if (.verbose) {
      message(sprintf(
        "Using dense SVD on scaled data with %i HVG.",
        length(selected_hvg)
      ))
    }
    zeallot::`%<-%`(
      c(pca_factors, pca_loadings, singular_values, scaled),
      rs_sc_pca(
        f_path_gene = get_rust_count_gene_f_path(object),
        f_path_cell = get_rust_count_cell_f_path(object),
        no_pcs = no_pcs,
        pca_params = pca_params,
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
  }

  if (.verbose) {
    message(sprintf(
      "Using sparse SVD on scaled data with %i HVG.",
      length(selected_hvg)
    ))
  }
  zeallot::`%<-%`(
    c(sparse_pca_factors, sparse_pca_loadings, sparse_pca_eigenvals),
    rs_sc_pca_sparse(
      f_path_gene = get_rust_count_gene_f_path(object),
      f_path_cell = get_rust_count_cell_f_path(object),
      no_pcs = no_pcs,
      pca_params = pca_params,
      cell_indices = get_cells_to_keep(object),
      gene_indices = selected_hvg,
      seed = seed,
      verbose = parse_verbosity(.verbose)
    )
  )

  object <- set_pca_factors(object, sparse_pca_factors)
  object <- set_pca_loadings(object, sparse_pca_loadings)
  object <- set_pca_singular_vals(object, sparse_pca_eigenvals[1:no_pcs])
  object
}

## fast clustering -------------------------------------------------------------

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

  checkmate::assertClass(object, "bixverse::SingleCellsSubset")
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

  # cell_idx is 1-indexed ORIGINAL positions; matches obs_table$cell_idx
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
