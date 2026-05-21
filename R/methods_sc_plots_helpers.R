# single cell plot helpers -----------------------------------------------------

## 2d embeddings ---------------------------------------------------------------

### helpers --------------------------------------------------------------------

#' Helper function to generate manifoldsR nearest neighbours
#'
#' @param x `SingleCells`, `MetaCells` object from which to extract the kNN
#' data.
#' @param modality String. One of `c("rna", "adt")`.
#'
#' @returns The manifoldsR `NearestNeigbours` if possible
#'
#' @keywords internal
.get_manifoldsr_knn <- function(x, modality = c("rna", "adt")) {
  modality <- match.arg(modality)
  checkmate::assertTRUE(
    S7::S7_inherits(x, SingleCells) || S7::S7_inherits(x, MetaCells)
  )

  knn_obj <- get_knn_obj(x, modality = modality)
  if (is.null(knn_obj)) {
    warning(paste(
      "No kNN data found.",
      "Attempting to use embeddings for 2D embedding"
    ))
    return(knn_obj)
  }
  manifold_nn <- sc_knn_to_nearest_neighbours(knn_obj)

  return(manifold_nn)
}

#' Placeholder: manifoldsR nearest neighbours from a WNN graph
#'
#' @keywords internal
.get_manifoldsr_knn_from_wnn <- function(x) {
  checkmate::assertTRUE(
    S7::S7_inherits(x, SingleCellsMultiModal)
  )

  res <- S7::prop(x, "other_data")[["wnn"]][["knn"]]

  if (is.null(res)) {
    stop("WNN-based neighbours were not found.")
  }

  manifold_nn <- sc_knn_to_nearest_neighbours(res)

  return(manifold_nn)
}

### umap -----------------------------------------------------------------------

#' Run UMAP on a SingleCells/MetaCells object
#'
#' @description
#' Wrapper around [manifoldsR::umap()] for the `SingleCells` and `MetaCells`
#' classes. UMAP produces a low-dimensional embedding that emphasises local
#' neighbourhood structure while being computationally efficient via its
#' negative-sampling-based optimisation. It is the de facto default for
#' visualising single-cell data, though claims that it preserves global
#' structure substantially better than t-SNE are not well supported; with
#' matched initialisation (e.g. PCA or Laplacian Eigenmaps), the two methods
#' behave similarly on global geometry, and both should be interpreted primarily
#' as views of local structure.
#'
#' When `use_knn = TRUE` (the default), the kNN graph already stored on the
#' object (via [bixverse::find_neighbours_sc()]) is reused, which avoids
#' recomputing nearest neighbours and keeps the UMAP consistent with any
#' downstream sNN-based clustering. If no kNN is present, neighbours are
#' computed from the chosen embedding on the fly.
#'
#' Key parameters to tune: `k` controls the balance between local and global
#' structure (larger values produce more global layouts), while `min_dist`
#' and `spread` together control how tightly points are packed in the
#' embedding. For `MetaCells`, smaller `k` values are often appropriate given
#' the reduced number of points.
#'
#' @param object `SingleCells`, `MetaCells` class.
#' @param use_knn Boolean. Use the kNN graph found in the object. Defaults to
#' `TRUE`. If not available, will default to the embedding.
#' @param embd_to_use String. The embedding to use for UMAP. Must be available
#' in the object.
#' @param slot_name String. The name of this embedding within the object.
#' Defaults to `"umap"`.
#' @param no_embd_to_use Optional integer. Number of embedding dimensions to
#' use. If `NULL` all will be used.
#' @param modality String. On which modality to run the UMAP. One of
#' `c("rna", "adt", "wnn")`. The two latter options are only available for
#' multi-modal versions with the added data.
#' @param n_dim Integer. Number of UMAP dimensions. Defaults to `2L`.
#' @param k Integer. Number of nearest neighbours. Defaults to `15L`.
#' @param min_dist Numeric. Minimum distance between embedded points. Defaults
#' to `0.5`.
#' @param spread Numeric. Effective scale of embedded points. Defaults to `1.0`.
#' @param knn_method String. Approximate nearest neighbour algorithm. One of
#' `"hnsw"`, `"balltree"`, `"annoy"`, `"nndescent"`, or `"exhaustive"`.
#' @param nn_params Named list. See [manifoldsR::params_nn()].
#' @param umap_params Named list. See [manifoldsR::params_umap()].
#' @param seed Integer. For reproducibility.
#' @param .verbose Boolean. Controls verbosity.
#'
#' @return The object with a `"umap"` embedding added.
#'
#' @export
umap_sc <- S7::new_generic(
  name = "umap_sc",
  dispatch_args = "object",
  fun = function(
    object,
    use_knn = TRUE,
    embd_to_use = "pca",
    slot_name = "umap",
    no_embd_to_use = NULL,
    modality = c("rna", "adt", "wnn"),
    n_dim = 2L,
    k = 15L,
    min_dist = 0.5,
    spread = 1.0,
    knn_method = c(
      "kmknn",
      "hnsw",
      "balltree",
      "annoy",
      "nndescent",
      "exhaustive"
    ),
    nn_params = manifoldsR::params_nn(),
    umap_params = manifoldsR::params_umap(),
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

S7::method(umap_sc, ScOrMc) <- function(
  object,
  use_knn = TRUE,
  embd_to_use = "pca",
  slot_name = "umap",
  no_embd_to_use = NULL,
  modality = c("rna", "adt", "wnn"),
  n_dim = 2L,
  k = 15L,
  min_dist = 0.5,
  spread = 1.0,
  knn_method = c(
    "kmknn",
    "hnsw",
    "balltree",
    "annoy",
    "nndescent",
    "exhaustive"
  ),
  nn_params = manifoldsR::params_nn(),
  umap_params = manifoldsR::params_umap(),
  seed = 42L,
  .verbose = TRUE
) {
  modality <- match.arg(modality)

  # checks
  checkmate::assertTRUE(
    S7::S7_inherits(object, SingleCells) || S7::S7_inherits(object, MetaCells)
  )
  checkmate::qassert(use_knn, "B1")
  checkmate::qassert(embd_to_use, "S1")
  checkmate::qassert(slot_name, "S1")
  checkmate::qassert(no_embd_to_use, c("I1", "0"))
  checkmate::qassert(n_dim, "I1[1,)")
  checkmate::qassert(k, "I1[2,)")
  checkmate::qassert(min_dist, "N1[0,)")
  checkmate::qassert(spread, "N1[0,)")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, "B1")
  knn_method <- match.arg(knn_method)

  if (modality != "rna" && !S7::S7_inherits(object, SingleCellsMultiModal)) {
    stop(sprintf(
      "modality = '%s' is only supported for SingleCellsMultiModal.",
      modality
    ))
  }

  # wnn takes the integrated graph; embeddings still read/write the rna cache
  cache_modality <- if (modality == "wnn") "rna" else modality

  # get the knn
  knn <- if (modality == "wnn") {
    .get_manifoldsr_knn_from_wnn(x = object)
  } else if (use_knn) {
    .get_manifoldsr_knn(x = object, modality = modality)
  } else {
    NULL
  }

  # get the embedding
  checkmate::assertTRUE(
    embd_to_use %in% get_available_embeddings(object, modality = cache_modality)
  )
  embd <- get_embedding(
    x = object,
    embd_name = embd_to_use,
    modality = cache_modality
  )

  if (!is.null(no_embd_to_use)) {
    to_take <- min(c(no_embd_to_use, ncol(embd)))
    embd <- embd[, 1:to_take]
  }

  if (.verbose) {
    message("Running UMAP.")
  }

  umap_embd <- manifoldsR::umap(
    data = embd,
    n_dim = n_dim,
    knn = knn,
    k = k,
    min_dist = min_dist,
    spread = spread,
    knn_method = knn_method,
    nn_params = nn_params,
    umap_params = umap_params,
    seed = seed,
    .verbose = .verbose
  )

  colnames(umap_embd) <- sprintf("umap_%s", seq_len(ncol(umap_embd)))

  slot <- ifelse(modality == "wnn", "other", cache_modality)

  object <- set_embedding(
    x = object,
    embd = umap_embd,
    name = slot_name,
    modality = slot
  )

  return(object)
}

### tsne -----------------------------------------------------------------------

#' Run t-SNE on a SingleCells/MetaCells object
#'
#' @description
#' Wrapper around [manifoldsR::tsne()] for the `SingleCells` and `MetaCells`
#' classes. t-SNE produces a low-dimensional embedding that emphasises local
#' neighbourhood structure. Distances between well-separated clusters should not
#' be over-interpreted quantitatively, but the common claim that t-SNE discards
#' global structure while UMAP preserves it is largely an artefact of default
#' initialisations rather than a property of the loss functions themselves.
#'
#' When `use_knn = FALSE` (the default), the kNN graph already stored on the
#' object is reused. Otherwise neighbours are computed from the chosen
#' embedding.
#'
#' Two approximation strategies are available via `approx_type`: `"bh"`
#' (Barnes-Hut) is the classical O(n log n) approximation and works well
#' across a wide range of dataset sizes; `"fft"` (interpolation-based, as in
#' FIt-SNE) scales better to very large datasets. `perplexity` controls the
#' bandwidth of the Gaussian kernel used to compute affinities within the
#' neighbour set (typical values 5-50). When a pre-computed kNN is supplied via
#' `use_knn = TRUE`, perplexity no longer drives neighbour retrieval but still
#' shapes the affinity distribution over the retrieved neighbours; values too
#' close to the kNN size will produce poor results. With tSNE in particular the
#' rule of thumb is to set k to `3 * perplexity`. When `k ≤ perplexity`` the
#' algorithm does not behave properly anymore, thus, will throw an error.
#'
#' @param object `SingleCells`, `MetaCells` class.
#' @param use_knn Boolean. Use the kNN graph found in the object. Defaults to
#' `TRUE`. If not available, will default to the embedding.
#' @param embd_to_use String. The embedding to use for t-SNE. Must be available
#' in the object.
#' @param slot_name String. The name of this embedding within the object.
#' Defaults to `"tsne"`.
#' @param no_embd_to_use Optional integer. Number of embedding dimensions to
#' use. If `NULL` all will be used.
#' @param modality String. On which modality to run the UMAP. One of
#' `c("rna", "adt", "wnn")`. The two latter options are only available for
#' multi-modal versions with the added data.
#' @param n_dim Integer. Number of t-SNE dimensions. Currently only `2L` is
#' supported. Defaults to `2L`.
#' @param perplexity Numeric. Perplexity parameter. Typical values between 5
#' and 50. Defaults to `30.0`.
#' @param approx_type String. Approximation method. One of `"bh"` (Barnes-Hut)
#' or `"fft"`. Defaults to `"bh"`.
#' @param knn_method String. Approximate nearest neighbour algorithm. One of
#' `"hnsw"`, `"balltree"`, `"annoy"`, `"nndescent"`, or `"exhaustive"`.
#' @param nn_params Named list. See [manifoldsR::params_nn()].
#' @param tsne_params Named list. See [manifoldsR::params_tsne()].
#' @param seed Integer. For reproducibility.
#' @param .verbose Boolean. Controls verbosity.
#'
#' @return The object with a `"tsne"` embedding added.
#'
#' @export
tsne_sc <- S7::new_generic(
  name = "tsne_sc",
  dispatch_args = "object",
  fun = function(
    object,
    use_knn = FALSE,
    embd_to_use = "pca",
    slot_name = "tsne",
    no_embd_to_use = NULL,
    modality = c("rna", "adt", "wnn"),
    n_dim = 2L,
    perplexity = 10.0,
    approx_type = c("bh", "fft"),
    knn_method = c(
      "kmknn",
      "hnsw",
      "balltree",
      "annoy",
      "nndescent",
      "exhaustive"
    ),
    nn_params = manifoldsR::params_nn(),
    tsne_params = manifoldsR::params_tsne(),
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

S7::method(tsne_sc, ScOrMc) <- function(
  object,
  use_knn = FALSE,
  embd_to_use = "pca",
  slot_name = "tsne",
  no_embd_to_use = NULL,
  modality = c("rna", "adt", "wnn"),
  n_dim = 2L,
  perplexity = 10.0,
  approx_type = c("bh", "fft"),
  knn_method = c(
    "kmknn",
    "hnsw",
    "balltree",
    "annoy",
    "nndescent",
    "exhaustive"
  ),
  nn_params = manifoldsR::params_nn(),
  tsne_params = manifoldsR::params_tsne(),
  seed = 42L,
  .verbose = TRUE
) {
  modality <- match.arg(modality)

  # checks
  checkmate::assertTRUE(
    S7::S7_inherits(object, SingleCells) || S7::S7_inherits(object, MetaCells)
  )
  checkmate::qassert(embd_to_use, "S1")
  checkmate::qassert(slot_name, "S1")
  checkmate::qassert(no_embd_to_use, c("I1", "0"))
  checkmate::qassert(n_dim, "I1[2,2]")
  checkmate::qassert(perplexity, "N1[1,)")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, "B1")
  approx_type <- match.arg(approx_type)
  knn_method <- match.arg(knn_method)

  if (modality != "rna" && !S7::S7_inherits(object, SingleCellsMultiModal)) {
    stop(sprintf(
      "modality = '%s' is only supported for SingleCellsMultiModal.",
      modality
    ))
  }

  cache_modality <- if (modality == "wnn") "rna" else modality

  # get the knn
  knn <- if (modality == "wnn") {
    .get_manifoldsr_knn_from_wnn(x = object)
  } else if (use_knn) {
    .get_manifoldsr_knn(x = object, modality = modality)
  } else {
    NULL
  }

  # get the embedding
  checkmate::assertTRUE(
    embd_to_use %in% get_available_embeddings(object, modality = cache_modality)
  )
  embd <- get_embedding(
    x = object,
    embd_name = embd_to_use,
    modality = cache_modality
  )

  if (!is.null(no_embd_to_use)) {
    to_take <- min(c(no_embd_to_use, ncol(embd)))
    embd <- embd[, 1:to_take]
  }

  if (.verbose) {
    message("Running t-SNE.")
  }

  tsne_embd <- manifoldsR::tsne(
    data = embd,
    knn = knn,
    n_dim = n_dim,
    perplexity = perplexity,
    approx_type = approx_type,
    knn_method = knn_method,
    nn_params = nn_params,
    tsne_params = tsne_params,
    seed = seed,
    .verbose = .verbose
  )

  colnames(tsne_embd) <- sprintf("tsne_%s", seq_len(ncol(tsne_embd)))

  slot <- ifelse(modality == "wnn", "other", cache_modality)

  object <- set_embedding(
    x = object,
    embd = tsne_embd,
    name = slot_name,
    modality = slot
  )

  return(object)
}

### phate ----------------------------------------------------------------------

#' Run PHATE on a SingleCells/MetaCells object
#'
#' @description
#' Wrapper around [manifoldsR::phate()] for the `SingleCells` and `MetaCells`
#' classes. PHATE (Potential of Heat-diffusion for Affinity-based Trajectory
#' Embedding) produces a low-dimensional embedding that preserves both local and
#' global structure by operating on a diffusion process over the data manifold.
#' Unlike UMAP or t-SNE, PHATE is explicitly designed to reveal continuous
#' progressions and branching structure, making it the preferred choice for data
#' with developmental or trajectory-like organisation.
#'
#' When `use_knn = TRUE` (the default), the kNN graph already stored on the
#' object is reused; otherwise neighbours are computed from the chosen
#' embedding. The algorithm then constructs a diffusion operator, raises it to a
#' power (the diffusion time `t`, see [manifoldsR::params_phate()]) that
#' denoises the manifold, and computes potential distances that are finally
#' embedded via metric MDS.
#'
#' Because PHATE inherently smooths over the kNN graph, it pairs naturally
#' with `MetaCells`: the combination yields a particularly clean view of
#' continuous biological processes on denoised data.
#'
#' @param object `SingleCells`, `MetaCells` class.
#' @param use_knn Boolean. Use the kNN graph found in the object. Defaults to
#' `TRUE`. If not available, will default to the embedding.
#' @param embd_to_use String. The embedding to use for PHATE. Must be available
#' in the object.
#' @param slot_name String. The name of this embedding within the object.
#' Defaults to `"phate"`.
#' @param no_embd_to_use Optional integer. Number of embedding dimensions to
#' use. If `NULL` all will be used.
#' @param modality String. On which modality to run the UMAP. One of
#' `c("rna", "adt", "wnn")`. The two latter options are only available for
#' multi-modal versions with the added data.
#' @param n_dim Integer. Number of PHATE dimensions. Currently only `2L` is
#' supported. Defaults to `2L`.
#' @param k Integer. Number of nearest neighbours for graph construction.
#' Defaults to `5L`.
#' @param knn_method String. Approximate nearest neighbour algorithm. One of
#' `"hnsw"`, `"balltree"`, `"annoy"`, `"nndescent"`, or `"exhaustive"`.
#' @param nn_params Named list. See [manifoldsR::params_nn()].
#' @param phate_params Named list. See [manifoldsR::params_phate()].
#' @param seed Integer. For reproducibility.
#' @param .verbose Boolean. Controls verbosity.
#'
#' @return The object with a `"phate"` embedding added.
#'
#' @export
phate_sc <- S7::new_generic(
  name = "phate_sc",
  dispatch_args = "object",
  fun = function(
    object,
    use_knn = TRUE,
    embd_to_use = "pca",
    slot_name = "phate",
    no_embd_to_use = NULL,
    modality = c("rna", "adt", "wnn"),
    n_dim = 2L,
    k = 5L,
    knn_method = c(
      "kmknn",
      "hnsw",
      "balltree",
      "annoy",
      "nndescent",
      "exhaustive"
    ),
    nn_params = manifoldsR::params_nn(),
    phate_params = manifoldsR::params_phate(),
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

S7::method(phate_sc, ScOrMc) <- function(
  object,
  use_knn = TRUE,
  embd_to_use = "pca",
  slot_name = "phate",
  no_embd_to_use = NULL,
  modality = c("rna", "adt", "wnn"),
  n_dim = 2L,
  k = 5L,
  knn_method = c(
    "kmknn",
    "hnsw",
    "balltree",
    "annoy",
    "nndescent",
    "exhaustive"
  ),
  nn_params = manifoldsR::params_nn(),
  phate_params = manifoldsR::params_phate(),
  seed = 42L,
  .verbose = TRUE
) {
  modality <- match.arg(modality)

  # checks
  checkmate::assertTRUE(
    S7::S7_inherits(object, SingleCells) || S7::S7_inherits(object, MetaCells)
  )
  checkmate::qassert(embd_to_use, "S1")
  checkmate::qassert(slot_name, "S1")
  checkmate::qassert(no_embd_to_use, c("I1", "0"))
  checkmate::qassert(n_dim, "I1[2,2]")
  checkmate::qassert(k, "I1[1,)")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, "B1")
  knn_method <- match.arg(knn_method)

  if (modality != "rna" && !S7::S7_inherits(object, SingleCellsMultiModal)) {
    stop(sprintf(
      "modality = '%s' is only supported for SingleCellsMultiModal.",
      modality
    ))
  }

  cache_modality <- if (modality == "wnn") "rna" else modality

  # get the knn
  knn <- if (modality == "wnn") {
    .get_manifoldsr_knn_from_wnn(x = object)
  } else if (use_knn) {
    .get_manifoldsr_knn(x = object, modality = modality)
  } else {
    NULL
  }

  # get the embedding
  checkmate::assertTRUE(
    embd_to_use %in% get_available_embeddings(object, modality = cache_modality)
  )
  embd <- get_embedding(
    x = object,
    embd_name = embd_to_use,
    modality = cache_modality
  )

  if (!is.null(no_embd_to_use)) {
    to_take <- min(c(no_embd_to_use, ncol(embd)))
    embd <- embd[, 1:to_take]
  }

  if (.verbose) {
    message("Running PHATE.")
  }

  phate_embd <- manifoldsR::phate(
    data = embd,
    knn = knn,
    n_dim = n_dim,
    k = k,
    knn_method = knn_method,
    nn_params = nn_params,
    phate_params = phate_params,
    seed = seed,
    .verbose = .verbose
  )

  colnames(phate_embd) <- sprintf("phate_%s", seq_len(ncol(phate_embd)))

  slot <- ifelse(modality == "wnn", "other", cache_modality)

  object <- set_embedding(
    x = object,
    embd = phate_embd,
    name = slot_name,
    modality = slot
  )

  return(object)
}

## feature extraction ----------------------------------------------------------

### helpers --------------------------------------------------------------------

#' Match requested features against available feature names
#'
#' @param features Character vector. Requested feature ids.
#' @param available Character vector. Feature names present in the data.
#'
#' @return The subset of `features` that was matched, in input order.
#'
#' @keywords internal
.match_features <- function(features, available) {
  idx <- match(features, available)
  missing <- is.na(idx)
  if (any(missing)) {
    warning(sprintf("%i features could not be matched.", sum(missing)))
    features <- features[!missing]
  }
  if (length(features) == 0) {
    stop("No features matched. Please double check provided parameters!")
  }
  features
}

#' Z-score and optionally clip a numeric vector
#'
#' @param x Numeric vector.
#' @param clip Optional numeric. Clip the z-scores to `[-clip, clip]`.
#'
#' @return The scaled (and optionally clipped) vector.
#'
#' @keywords internal
.scale_and_clip_expr <- function(x, clip = NULL) {
  s <- stats::sd(x)
  if (s > 1e-8) {
    x <- (x - mean(x)) / s
  }
  if (!is.null(clip)) {
    x <- pmax(pmin(x, clip), -clip)
  }
  x
}

#' Extract dense expression for in-memory count matrices
#'
#' @param mat Numeric matrix. Cells x features, with cell ids as row names and
#' feature ids as column names.
#' @param features Character vector. Feature ids to extract (pre-matched).
#' @param scale Boolean. Whether to z-score each feature across cells.
#' @param clip Optional numeric. Clip the z-scores if `scale = TRUE`.
#'
#' @return A data.table with a `cell_id` column and one column per feature.
#'
#' @keywords internal
.extract_expr_in_memory <- function(mat, features, scale = FALSE, clip = NULL) {
  idx <- match(features, colnames(mat))
  dt <- data.table::data.table(cell_id = rownames(mat))
  for (i in seq_along(features)) {
    v <- mat[, idx[i]]
    if (scale) {
      v <- .scale_and_clip_expr(v, clip)
    }
    data.table::set(dt, j = features[i], value = v)
  }
  dt
}

#' Per-group mean expression and percent expressed for in-memory matrices
#'
#' @param mat Numeric matrix. Cells x features.
#' @param features Character vector. Feature ids to extract (pre-matched).
#' @param grouping Factor. Group assignment per cell.
#'
#' @return A long data.table with columns `gene`, `group`, `mean_exp`,
#' `pct_exp`.
#'
#' @keywords internal
.grouped_gene_stats_in_memory <- function(mat, features, grouping) {
  group_levels <- levels(grouping)
  sub <- mat[, match(features, colnames(mat)), drop = FALSE]

  res <- lapply(seq_along(features), function(j) {
    vals <- sub[, j]
    data.table::data.table(
      gene = features[j],
      group = group_levels,
      mean_exp = as.numeric(tapply(vals, grouping, mean)),
      pct_exp = as.numeric(tapply(vals > 0, grouping, mean)) * 100
    )
  })

  data.table::rbindlist(res)
}

#' Finalise a long dot-plot data.table
#'
#' @param plot_dt data.table. With columns `gene`, `group`, `mean_exp`,
#' `pct_exp`.
#' @param features Character vector. Feature ids, in display order.
#' @param group_levels Character vector. Group levels, in display order.
#' @param scale_exp Boolean. Whether to min-max scale mean expression per gene.
#'
#' @return The data.table with an added `scaled_exp` column and ordered
#' `gene`/`group` factors.
#'
#' @keywords internal
.finalise_dot_plot_dt <- function(plot_dt, features, group_levels, scale_exp) {
  plot_dt[, scaled_exp := mean_exp]
  if (scale_exp) {
    plot_dt[,
      scaled_exp := {
        rng <- range(mean_exp)
        if (rng[1] == rng[2]) 0 else (mean_exp - rng[1]) / (rng[2] - rng[1])
      },
      by = gene
    ]
  }
  plot_dt[, gene := factor(gene, levels = features)]
  plot_dt[, group := factor(group, levels = group_levels)]
  plot_dt[]
}

### gene summaries -------------------------------------------------------------

# generic in base_generics_sc.R

#' @method extract_dot_plot_data SingleCells
#'
#' @export
S7::method(extract_dot_plot_data, SingleCells) <- function(
  object,
  features,
  grouping_variable,
  scale_exp = TRUE,
  modality = c("rna", "adt")
) {
  modality <- match.arg(modality)
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::qassert(features, "S+")
  checkmate::qassert(grouping_variable, "S1")
  checkmate::qassert(scale_exp, "B1")

  if (modality != "rna") {
    stop(paste(
      "SingleCells only supports modality = 'rna'.",
      "Use SingleCellsMultiModal for ADT."
    ))
  }

  gene_idx <- get_gene_indices(
    x = object,
    gene_ids = features,
    rust_index = TRUE
  )

  grouping <- as.factor(
    unlist(object[[grouping_variable]], use.names = FALSE)
  )

  gene_res <- rs_extract_grouped_gene_stats(
    f_path = get_rust_count_gene_f_path(object),
    cell_indices = get_cells_to_keep(object),
    gene_indices = gene_idx,
    group_ids = as.integer(grouping) - 1L,
    group_levels = levels(grouping)
  )

  n_groups <- length(gene_res$grp_label)

  plot_dt <- data.table::data.table(
    gene = rep(features, each = n_groups),
    group = rep(gene_res$grp_label, times = length(features)),
    mean_exp = gene_res$mean_exp,
    pct_exp = gene_res$perc_exp * 100
  )

  .finalise_dot_plot_dt(plot_dt, features, levels(grouping), scale_exp)
}

#' @method extract_dot_plot_data SingleCellsMultiModal
#'
#' @export
S7::method(extract_dot_plot_data, SingleCellsMultiModal) <- function(
  object,
  features,
  grouping_variable,
  scale_exp = TRUE,
  modality = c("rna", "adt")
) {
  modality <- match.arg(modality)

  if (modality == "rna") {
    rna_method <- S7::method(extract_dot_plot_data, SingleCells)
    return(rna_method(
      object = object,
      features = features,
      grouping_variable = grouping_variable,
      scale_exp = scale_exp,
      modality = "rna"
    ))
  }

  # ADT path
  checkmate::qassert(features, "S+")
  checkmate::qassert(grouping_variable, "S1")
  checkmate::qassert(scale_exp, "B1")

  adt <- S7::prop(object, "adt_counts")
  if (is.null(adt)) {
    stop("No ADT counts in this object. Add them with add_adt_counts_sc().")
  }

  features <- .match_features(features, colnames(adt$norm_counts))
  grouping <- as.factor(
    unlist(object[[grouping_variable]], use.names = FALSE)
  )

  plot_dt <- .grouped_gene_stats_in_memory(adt$norm_counts, features, grouping)
  .finalise_dot_plot_dt(plot_dt, features, levels(grouping), scale_exp)
}

#' @method extract_dot_plot_data MetaCells
#'
#' @export
S7::method(extract_dot_plot_data, MetaCells) <- function(
  object,
  features,
  grouping_variable,
  scale_exp = TRUE,
  modality = c("rna", "adt")
) {
  modality <- match.arg(modality)
  checkmate::assertTRUE(S7::S7_inherits(object, MetaCells))
  checkmate::qassert(features, "S+")
  checkmate::qassert(grouping_variable, "S1")
  checkmate::qassert(scale_exp, "B1")

  if (modality != "rna") {
    stop(paste(
      "MetaCells only supports modality = 'rna'.",
      "Use SingleCellsMultiModal for ADT."
    ))
  }

  mat <- get_sc_counts(object, assay = "norm")
  features <- .match_features(features, colnames(mat))
  grouping <- as.factor(
    unlist(object[[grouping_variable]], use.names = FALSE)
  )

  plot_dt <- .grouped_gene_stats_in_memory(mat, features, grouping)
  .finalise_dot_plot_dt(plot_dt, features, levels(grouping), scale_exp)
}

### individual cells -----------------------------------------------------------

#' @method extract_gene_expression SingleCells
#'
#' @export
S7::method(extract_gene_expression, SingleCells) <- function(
  object,
  features,
  obs_cols = NULL,
  scale = FALSE,
  clip = NULL,
  modality = c("rna", "adt")
) {
  modality <- match.arg(modality)
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::qassert(features, "S+")
  checkmate::qassert(obs_cols, c("0", "S+"))
  checkmate::qassert(scale, "B1")
  checkmate::qassert(clip, c("0", "N1(0,)"))

  if (modality != "rna") {
    stop(paste(
      "SingleCells only supports modality = 'rna'.",
      "Use SingleCellsMultiModal for ADT."
    ))
  }

  gene_idx <- get_gene_indices(
    x = object,
    gene_ids = features,
    rust_index = TRUE
  )

  counts <- rs_extract_several_genes_plots(
    f_path = get_rust_count_gene_f_path(object),
    cell_indices = get_cells_to_keep(object),
    gene_indices = gene_idx,
    scale = scale,
    clip = clip
  )

  dt <- data.table::data.table(
    cell_id = get_cell_names(object, filtered = TRUE)
  )
  for (i in seq_along(features)) {
    data.table::set(dt, j = features[i], value = counts[[i]])
  }

  if (!is.null(obs_cols)) {
    obs_dt <- object[[obs_cols]]
    for (col in names(obs_dt)) {
      data.table::set(dt, j = col, value = obs_dt[[col]])
    }
  }

  dt
}

#' @method extract_gene_expression SingleCellsMultiModal
#'
#' @export
S7::method(extract_gene_expression, SingleCellsMultiModal) <- function(
  object,
  features,
  obs_cols = NULL,
  scale = FALSE,
  clip = NULL,
  modality = c("rna", "adt")
) {
  modality <- match.arg(modality)

  if (modality == "rna") {
    rna_method <- S7::method(extract_gene_expression, SingleCells)
    return(rna_method(
      object = object,
      features = features,
      obs_cols = obs_cols,
      scale = scale,
      clip = clip,
      modality = "rna"
    ))
  }

  # ADT path
  checkmate::qassert(features, "S+")
  checkmate::qassert(obs_cols, c("0", "S+"))
  checkmate::qassert(scale, "B1")
  checkmate::qassert(clip, c("0", "N1(0,)"))

  adt <- S7::prop(object, "adt_counts")
  if (is.null(adt)) {
    stop("No ADT counts in this object. Add them with add_adt_counts_sc().")
  }

  features <- .match_features(features, colnames(adt$norm_counts))
  dt <- .extract_expr_in_memory(adt$norm_counts, features, scale, clip)

  # obs are attached positionally; rows are the kept cells in kept order
  if (!is.null(obs_cols)) {
    obs_dt <- object[[obs_cols]]
    for (col in names(obs_dt)) {
      data.table::set(dt, j = col, value = obs_dt[[col]])
    }
  }

  dt
}

#' @method extract_gene_expression MetaCells
#'
#' @export
S7::method(extract_gene_expression, MetaCells) <- function(
  object,
  features,
  obs_cols = NULL,
  scale = FALSE,
  clip = NULL,
  modality = c("rna", "adt")
) {
  modality <- match.arg(modality)
  checkmate::assertTRUE(S7::S7_inherits(object, MetaCells))
  checkmate::qassert(features, "S+")
  checkmate::qassert(obs_cols, c("0", "S+"))
  checkmate::qassert(scale, "B1")
  checkmate::qassert(clip, c("0", "N1(0,)"))

  if (modality != "rna") {
    stop(paste(
      "MetaCells only supports modality = 'rna'.",
      "Use SingleCellsMultiModal for ADT."
    ))
  }

  mat <- get_sc_counts(object, assay = "norm")
  features <- .match_features(features, colnames(mat))
  dt <- .extract_expr_in_memory(mat, features, scale, clip)

  if (!is.null(obs_cols)) {
    obs_dt <- object[[obs_cols]]
    for (col in names(obs_dt)) {
      data.table::set(dt, j = col, value = obs_dt[[col]])
    }
  }

  dt
}
