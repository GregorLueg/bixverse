# single cell plot helpers -----------------------------------------------------

## 2d embeddings ---------------------------------------------------------------

### helpers --------------------------------------------------------------------

#' Helper function to generate manifoldsR nearest neighbours
#'
#' @param x `SingleCells`, `MetaCells` object from which to extract the kNN
#' data.
#'
#' @returns The manifoldsR `NearestNeigbours` if possible
#'
#' @keywords internal
.get_manifoldsr_knn <- function(x) {
  # checks
  checkmate::assertTRUE(
    S7::S7_inherits(x, SingleCells) || S7::S7_inherits(x, MetaCells)
  )

  knn_obj <- get_knn_obj(x)
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

### umap -----------------------------------------------------------------------

#' Run UMAP on a SingleCells object
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
#' @param no_embd_to_use Optional integer. Number of embedding dimensions to
#' use. If `NULL` all will be used.
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
    no_embd_to_use = NULL,
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
  no_embd_to_use = NULL,
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
  # checks
  checkmate::assertTRUE(
    S7::S7_inherits(object, SingleCells) || S7::S7_inherits(object, MetaCells)
  )
  checkmate::qassert(use_knn, "B1")
  checkmate::qassert(embd_to_use, "S1")
  checkmate::qassert(no_embd_to_use, c("I1", "0"))
  checkmate::qassert(n_dim, "I1[1,)")
  checkmate::qassert(k, "I1[2,)")
  checkmate::qassert(min_dist, "N1[0,)")
  checkmate::qassert(spread, "N1[0,)")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, "B1")
  knn_method <- match.arg(knn_method)

  # get the knn
  knn <- if (use_knn) {
    .get_manifoldsr_knn(x = object)
  } else {
    NULL
  }

  # get the embedding
  checkmate::assertTRUE(embd_to_use %in% get_available_embeddings(object))
  embd <- get_embedding(x = object, embd_name = embd_to_use)

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
  object <- set_embedding(x = object, embd = umap_embd, name = "umap")

  return(object)
}

### tsne -----------------------------------------------------------------------

#' Run t-SNE on a SingleCells object
#'
#' @description
#' Wrapper around [manifoldsR::tsne()] for the `SingleCells` and `MetaCells`
#' classes. t-SNE produces a low-dimensional embedding that emphasises local
#' neighbourhood structure. Distances between well-separated clusters should not
#' be over-interpreted quantitatively, but the common claim that t-SNE discards
#' global structure while UMAP preserves it is largely an artefact of default
#' initialisations rather than a property of the loss functions themselves.
#'
#' When `use_knn = TRUE` (the default), the kNN graph already stored on the
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
#' close to the kNN size will produce poor results.
#'
#' @param object `SingleCells`, `MetaCells` class.
#' @param use_knn Boolean. Use the kNN graph found in the object. Defaults to
#' `TRUE`. If not available, will default to the embedding.
#' @param embd_to_use String. The embedding to use for t-SNE. Must be available
#' in the object.
#' @param no_embd_to_use Optional integer. Number of embedding dimensions to
#' use. If `NULL` all will be used.
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
    use_knn = TRUE,
    embd_to_use = "pca",
    no_embd_to_use = NULL,
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
  use_knn = TRUE,
  embd_to_use = "pca",
  no_embd_to_use = NULL,
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
  # checks
  checkmate::assertTRUE(
    S7::S7_inherits(object, SingleCells) || S7::S7_inherits(object, MetaCells)
  )
  checkmate::qassert(embd_to_use, "S1")
  checkmate::qassert(no_embd_to_use, c("I1", "0"))
  checkmate::qassert(n_dim, "I1[2,2]")
  checkmate::qassert(perplexity, "N1[1,)")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, "B1")
  approx_type <- match.arg(approx_type)
  knn_method <- match.arg(knn_method)

  # get the knn
  knn <- if (use_knn) {
    .get_manifoldsr_knn(x = object)
  } else {
    NULL
  }

  # get the embedding
  checkmate::assertTRUE(embd_to_use %in% get_available_embeddings(object))
  embd <- get_embedding(x = object, embd_name = embd_to_use)

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
  object <- set_embedding(x = object, embd = tsne_embd, name = "tsne")

  return(object)
}

### phate ----------------------------------------------------------------------

#' Run PHATE on a SingleCells object
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
#' @param no_embd_to_use Optional integer. Number of embedding dimensions to
#' use. If `NULL` all will be used.
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
    no_embd_to_use = NULL,
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
  no_embd_to_use = NULL,
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
  # checks
  checkmate::assertTRUE(
    S7::S7_inherits(object, SingleCells) || S7::S7_inherits(object, MetaCells)
  )
  checkmate::qassert(embd_to_use, "S1")
  checkmate::qassert(no_embd_to_use, c("I1", "0"))
  checkmate::qassert(n_dim, "I1[2,2]")
  checkmate::qassert(k, "I1[1,)")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, "B1")
  knn_method <- match.arg(knn_method)

  # get the knn
  knn <- if (use_knn) {
    .get_manifoldsr_knn(x = object)
  } else {
    NULL
  }

  # get the embedding
  checkmate::assertTRUE(embd_to_use %in% get_available_embeddings(object))
  embd <- get_embedding(x = object, embd_name = embd_to_use)

  if (!is.null(no_embd_to_use)) {
    to_take <- min(c(no_embd_to_use, ncol(embd)))
    embd <- embd[, 1:to_take]
  }

  if (.verbose) {
    message("Running PHATE.")
  }

  phate_embd <- manifoldsR::phate(
    data = embd,
    n_dim = n_dim,
    k = k,
    knn_method = knn_method,
    nn_params = nn_params,
    phate_params = phate_params,
    seed = seed,
    .verbose = .verbose
  )

  colnames(phate_embd) <- sprintf("phate_%s", seq_len(ncol(phate_embd)))
  object <- set_embedding(x = object, embd = phate_embd, name = "phate")

  return(object)
}

## gene extracters -------------------------------------------------------------

### gene summaries -------------------------------------------------------------

#' Extract grouped gene statistics for dot plots
#'
#' @description
#' Extracts per-group mean expression and percentage of expressing cells for a
#' set of genes. Returns a long-format data.table suitable for dot plots.
#'
#' @param object `SingleCells` class.
#' @param features Character vector. Gene IDs to extract.
#' @param grouping_variable String. Column name in the obs table to group by.
#' @param scale_exp Boolean. Whether to min-max scale mean expression per gene.
#'
#' @return A data.table with columns: gene, group, mean_exp, scaled_exp, pct_exp.
#'
#' @export
extract_dot_plot_data <- S7::new_generic(
  name = "extract_dot_plot_data",
  dispatch_args = "object",
  fun = function(object, features, grouping_variable, scale_exp = TRUE) {
    S7::S7_dispatch()
  }
)

#' @method extract_dot_plot_data SingleCells
#'
#' @export
S7::method(extract_dot_plot_data, SingleCells) <- function(
  object,
  features,
  grouping_variable,
  scale_exp = TRUE
) {
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::qassert(features, "S+")
  checkmate::qassert(grouping_variable, "S1")
  checkmate::qassert(scale_exp, "B1")

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

  plot_dt[, gene := factor(gene, levels = rev(features))]
  plot_dt[, group := factor(group, levels = levels(grouping))]

  plot_dt
}

### individual cells -----------------------------------------------------------

#' Extract normalised gene expression for plotting
#'
#' @description
#' Extracts dense normalised (log1p) expression values for a set of genes,
#' optionally with additional observation metadata columns.
#'
#' @param object `SingleCells` class.
#' @param features Character vector. Gene IDs to extract.
#' @param obs_cols Optional character vector. Column names from the obs table
#' to include.
#' @param scale Boolean. Whether to z-score the expression values.
#' @param clip Optional numeric. If `scale = TRUE`, clip z-scores to
#' `[-clip, clip]`.
#'
#' @return A data.table with a `cell_id` column, one column per gene, and
#'   any requested obs columns.
#'
#' @export
extract_gene_expression <- S7::new_generic(
  name = "extract_gene_expression",
  dispatch_args = "object",
  fun = function(
    object,
    features,
    obs_cols = NULL,
    scale = FALSE,
    clip = NULL
  ) {
    S7::S7_dispatch()
  }
)

#' @method extract_gene_expression SingleCells
#'
#' @export
S7::method(extract_gene_expression, SingleCells) <- function(
  object,
  features,
  obs_cols = NULL,
  scale = FALSE,
  clip = NULL
) {
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::qassert(features, "S+")
  checkmate::qassert(obs_cols, c("0", "S+"))
  checkmate::qassert(scale, "B1")
  checkmate::qassert(clip, c("0", "N1(0,)"))

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
