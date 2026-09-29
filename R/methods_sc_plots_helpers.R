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
    S7::S7_inherits(x, SingleCells) ||
      S7::S7_inherits(x, MetaCells) ||
      S7::S7_inherits(x, SingleCellsSubset)
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
#' @param x `SingleCellsMultiModal` object holding a WNN graph.
#'
#' @returns A `manifoldsR` nearest neighbours object.
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

#' Resolve the kNN and input embedding for a manifold method
#'
#' @description
#' Shared entry point of the `*_sc` 2D embedding methods. `"wnn"` takes the
#' integrated kNN graph but reads the input embedding from the RNA cache, so the
#' modality the embedding comes from (`cache_modality`) can differ from the one
#' the result is written to.
#'
#' @param object `SingleCells`, `MetaCells` or `SingleCellsSubset` class.
#' @param use_knn Boolean. Use the kNN graph found in the object.
#' @param embd_to_use String. The embedding to feed the method. Must be
#' available in the object.
#' @param no_embd_to_use Optional integer. Number of embedding dimensions to
#' use. If `NULL` all will be used.
#' @param modality String. One of `c("rna", "adt", "wnn")`.
#'
#' @returns A list with
#' \itemize{
#'   \item knn - The manifoldsR `NearestNeighbours`, or `NULL` if the method
#'   should build its own.
#'   \item embd - The input embedding, cells x dimensions.
#'   \item cache_modality - The modality the embedding was read from.
#' }
#'
#' @keywords internal
.manifold_inputs <- function(
  object,
  use_knn,
  embd_to_use,
  no_embd_to_use,
  modality
) {
  # checks
  checkmate::assertTRUE(
    S7::S7_inherits(object, SingleCells) ||
      S7::S7_inherits(object, MetaCells) ||
      S7::S7_inherits(object, SingleCellsSubset)
  )
  checkmate::qassert(use_knn, "B1")
  checkmate::qassert(embd_to_use, "S1")
  checkmate::qassert(no_embd_to_use, c("I1", "0"))
  checkmate::assertChoice(modality, c("rna", "adt", "wnn"))

  if (modality != "rna" && !S7::S7_inherits(object, SingleCellsMultiModal)) {
    stop(sprintf(
      "modality = '%s' is only supported for SingleCellsMultiModal.",
      modality
    ))
  }

  # wnn takes the integrated graph; embeddings still read/write the rna cache
  cache_modality <- if (modality == "wnn") "rna" else modality

  # hard tier: the manifold is written back onto the object, and it is read
  # from `cache_modality` while the kNN comes from `modality`
  assert_sc_state(object, artefacts = embd_to_use, modality = cache_modality)
  if (modality == "wnn" || use_knn) {
    assert_sc_state(object, artefacts = "knn", modality = modality)
  }

  knn <- if (modality == "wnn") {
    .get_manifoldsr_knn_from_wnn(x = object)
  } else if (use_knn) {
    .get_manifoldsr_knn(x = object, modality = modality)
  } else {
    NULL
  }

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

  list(knn = knn, embd = embd, cache_modality = cache_modality)
}

#' Name a manifold embedding and write it back onto the object
#'
#' @param object `SingleCells`, `MetaCells` or `SingleCellsSubset` class.
#' @param embd Numeric matrix. The embedding, cells x dimensions.
#' @param prefix String. Column name prefix, e.g. `"umap"`.
#' @param slot_name String. Name of the embedding within the object.
#' @param modality String. Modality the embedding is written to.
#' @param from Character vector. Parent artefact names for the provenance
#' stamp, see [.manifold_from()].
#'
#' @returns The object with the embedding added.
#'
#' @keywords internal
.store_manifold <- function(object, embd, prefix, slot_name, modality, from) {
  checkmate::assertMatrix(embd, mode = "numeric")
  checkmate::qassert(prefix, "S1")
  checkmate::qassert(slot_name, "S1")
  checkmate::qassert(modality, "S1")
  checkmate::qassert(from, "S+")

  colnames(embd) <- sprintf("%s_%s", prefix, seq_len(ncol(embd)))

  set_embedding(
    x = object,
    embd = embd,
    name = slot_name,
    modality = modality,
    from = from
  )
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
#' @returns The object with a `"umap"` embedding added.
#'
#' @export
#'
#' @examples
#' # UMAP off the cached kNN graph
#' sc <- demo_single_cells()
#' sc <- umap_sc(sc, .verbose = FALSE)
#' dim(get_embedding(sc, "umap"))
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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
  knn_method <- match.arg(knn_method)

  # checks
  checkmate::qassert(slot_name, "S1")
  checkmate::qassert(n_dim, "I1[1,)")
  checkmate::qassert(k, "I1[2,)")
  checkmate::qassert(min_dist, "N1[0,)")
  checkmate::qassert(spread, "N1[0,)")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, "B1")

  inputs <- .manifold_inputs(
    object = object,
    use_knn = use_knn,
    embd_to_use = embd_to_use,
    no_embd_to_use = no_embd_to_use,
    modality = modality
  )

  if (.verbose) {
    message("Running UMAP.")
  }

  umap_embd <- manifoldsR::umap(
    data = inputs$embd,
    n_dim = n_dim,
    knn = inputs$knn,
    k = k,
    min_dist = min_dist,
    spread = spread,
    knn_method = knn_method,
    nn_params = nn_params,
    umap_params = umap_params,
    seed = seed,
    .verbose = .verbose
  )

  .store_manifold(
    object = object,
    embd = umap_embd,
    prefix = "umap",
    slot_name = slot_name,
    modality = modality,
    from = .manifold_from(
      embd_to_use = embd_to_use,
      cache_modality = inputs$cache_modality,
      modality = modality,
      has_knn = !is.null(inputs$knn)
    )
  )
}

### densmap --------------------------------------------------------------------

#' Run densMAP on a SingleCells/MetaCells object
#'
#' @description
#' Wrapper around [manifoldsR::densmap()] for the `SingleCells` and `MetaCells`
#' classes. densMAP is UMAP plus a density-preserving term: a tight population
#' stays tight in the embedding and a diffuse one stays diffuse. With plain
#' UMAP the relative size of a cluster on the plot tells you nothing, with
#' densMAP it does. Setting `lambda = 0` in [manifoldsR::params_densmap()]
#' gives you back plain UMAP.
#'
#' Neighbour handling is the same as in [bixverse::umap_sc()]: with
#' `use_knn = TRUE` (the default) the cached kNN graph is reused, otherwise
#' neighbours are computed from the chosen embedding.
#'
#' @param object `SingleCells`, `MetaCells` class.
#' @param use_knn Boolean. Use the kNN graph found in the object. Defaults to
#' `TRUE`. If not available, will default to the embedding.
#' @param embd_to_use String. The embedding to use for densMAP. Must be
#' available in the object.
#' @param slot_name String. The name of this embedding within the object.
#' Defaults to `"densmap"`.
#' @param no_embd_to_use Optional integer. Number of embedding dimensions to
#' use. If `NULL` all will be used.
#' @param modality String. On which modality to run densMAP. One of
#' `c("rna", "adt", "wnn")`. The two latter options are only available for
#' multi-modal versions with the added data.
#' @param n_dim Integer. Number of densMAP dimensions. Defaults to `2L`.
#' @param k Integer. Number of nearest neighbours. Defaults to `15L`.
#' @param min_dist Numeric. Minimum distance between embedded points. Defaults
#' to `0.5`.
#' @param spread Numeric. Effective scale of embedded points. Defaults to `1.0`.
#' @param knn_method String. Approximate nearest neighbour algorithm. One of
#' `"hnsw"`, `"balltree"`, `"annoy"`, `"nndescent"`, or `"exhaustive"`.
#' @param nn_params Named list. See [manifoldsR::params_nn()].
#' @param umap_params Named list. See [manifoldsR::params_umap()].
#' @param dens_params Named list. The density knobs, see
#' [manifoldsR::params_densmap()].
#' @param seed Integer. For reproducibility.
#' @param .verbose Boolean. Controls verbosity.
#'
#' @returns The object with a `"densmap"` embedding added.
#'
#' @export
#'
#' @references Narayan, Berger & Cho, Nat. Biotechnol., 2021
#'
#' @examples
#' # densMAP off the cached kNN graph
#' sc <- demo_single_cells()
#' sc <- densmap_sc(sc, .verbose = FALSE)
#' dim(get_embedding(sc, "densmap"))
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
densmap_sc <- S7::new_generic(
  name = "densmap_sc",
  dispatch_args = "object",
  fun = function(
    object,
    use_knn = TRUE,
    embd_to_use = "pca",
    slot_name = "densmap",
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
    dens_params = manifoldsR::params_densmap(),
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

S7::method(densmap_sc, ScOrMc) <- function(
  object,
  use_knn = TRUE,
  embd_to_use = "pca",
  slot_name = "densmap",
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
  dens_params = manifoldsR::params_densmap(),
  seed = 42L,
  .verbose = TRUE
) {
  modality <- match.arg(modality)
  knn_method <- match.arg(knn_method)

  # checks
  checkmate::qassert(slot_name, "S1")
  checkmate::qassert(n_dim, "I1[1,)")
  checkmate::qassert(k, "I1[2,)")
  checkmate::qassert(min_dist, "N1[0,)")
  checkmate::qassert(spread, "N1[0,)")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, "B1")

  inputs <- .manifold_inputs(
    object = object,
    use_knn = use_knn,
    embd_to_use = embd_to_use,
    no_embd_to_use = no_embd_to_use,
    modality = modality
  )

  if (.verbose) {
    message("Running densMAP.")
  }

  densmap_embd <- manifoldsR::densmap(
    data = inputs$embd,
    knn = inputs$knn,
    n_dim = n_dim,
    k = k,
    min_dist = min_dist,
    spread = spread,
    knn_method = knn_method,
    nn_params = nn_params,
    umap_params = umap_params,
    dens_params = dens_params,
    seed = seed,
    .verbose = .verbose
  )

  .store_manifold(
    object = object,
    embd = densmap_embd,
    prefix = "densmap",
    slot_name = slot_name,
    modality = modality,
    from = .manifold_from(
      embd_to_use = embd_to_use,
      cache_modality = inputs$cache_modality,
      modality = modality,
      has_knn = !is.null(inputs$knn)
    )
  )
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
#' When `use_knn = TRUE`, the kNN graph already stored on the object is reused.
#' With `use_knn = FALSE` (the default), neighbours are computed from the chosen
#' embedding.
#'
#' Three approximation strategies are available via `approx_type`: `"bh"`
#' (Barnes-Hut) is the classical O(n log n) approximation and works well
#' across a wide range of dataset sizes; `"fft"` (interpolation-based, as in
#' FIt-SNE) scales better to very large datasets; `"fft_3k"` is the
#' three-kernel variant of the latter, with one forward and three inverse FFTs
#' per epoch instead of four each. The FFT options are only available on Unix
#' systems. `perplexity` controls the bandwidth of the Gaussian kernel used to
#' compute affinities within the neighbour set (typical values 5-50). When a
#' pre-computed kNN is supplied via `use_knn = TRUE`, perplexity no longer
#' drives neighbour retrieval but still shapes the affinity distribution over
#' the retrieved neighbours; values too close to the kNN size will produce poor
#' results. With tSNE in particular the rule of thumb is to set k to
#' `3 * perplexity`. When `k <= perplexity` the algorithm does not behave
#' properly anymore, thus, will throw an error.
#'
#' @param object `SingleCells`, `MetaCells` class.
#' @param use_knn Boolean. Use the kNN graph found in the object. Defaults to
#' `FALSE`. If not available, will default to the embedding.
#' @param embd_to_use String. The embedding to use for t-SNE. Must be available
#' in the object.
#' @param slot_name String. The name of this embedding within the object.
#' Defaults to `"tsne"`.
#' @param no_embd_to_use Optional integer. Number of embedding dimensions to
#' use. If `NULL` all will be used.
#' @param modality String. On which modality to run the t-SNE. One of
#' `c("rna", "adt", "wnn")`. The two latter options are only available for
#' multi-modal versions with the added data.
#' @param n_dim Integer. Number of t-SNE dimensions. Currently only `2L` is
#' supported. Defaults to `2L`.
#' @param perplexity Numeric. Perplexity parameter. Typical values between 5
#' and 50. Defaults to `10.0`.
#' @param approx_type String. Approximation method. One of `"bh"` (Barnes-Hut),
#' `"fft"` or `"fft_3k"`. Defaults to `"bh"`.
#' @param knn_method String. Approximate nearest neighbour algorithm. One of
#' `"hnsw"`, `"balltree"`, `"annoy"`, `"nndescent"`, or `"exhaustive"`.
#' @param nn_params Named list. See [manifoldsR::params_nn()].
#' @param tsne_params Named list. See [manifoldsR::params_tsne()].
#' @param seed Integer. For reproducibility.
#' @param .verbose Boolean. Controls verbosity.
#'
#' @returns The object with a `"tsne"` embedding added.
#'
#' @export
#'
#' @examples
#' # Barnes-Hut t-SNE on the PCA factors
#' sc <- demo_single_cells()
#' sc <- tsne_sc(sc, .verbose = FALSE)
#' dim(get_embedding(sc, "tsne"))
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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
    approx_type = c("bh", "fft", "fft_3k"),
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
  approx_type = c("bh", "fft", "fft_3k"),
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
  approx_type <- match.arg(approx_type)
  knn_method <- match.arg(knn_method)

  # checks
  checkmate::qassert(slot_name, "S1")
  checkmate::qassert(n_dim, "I1[2,2]")
  checkmate::qassert(perplexity, "N1[1,)")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, "B1")

  inputs <- .manifold_inputs(
    object = object,
    use_knn = use_knn,
    embd_to_use = embd_to_use,
    no_embd_to_use = no_embd_to_use,
    modality = modality
  )

  if (.verbose) {
    message("Running t-SNE.")
  }

  tsne_embd <- manifoldsR::tsne(
    data = inputs$embd,
    knn = inputs$knn,
    n_dim = n_dim,
    perplexity = perplexity,
    approx_type = approx_type,
    knn_method = knn_method,
    nn_params = nn_params,
    tsne_params = tsne_params,
    seed = seed,
    .verbose = .verbose
  )

  .store_manifold(
    object = object,
    embd = tsne_embd,
    prefix = "tsne",
    slot_name = slot_name,
    modality = modality,
    from = .manifold_from(
      embd_to_use = embd_to_use,
      cache_modality = inputs$cache_modality,
      modality = modality,
      has_knn = !is.null(inputs$knn)
    )
  )
}

### densne ---------------------------------------------------------------------

#' Run den-SNE on a SingleCells/MetaCells object
#'
#' @description
#' Wrapper around [manifoldsR::densne()] for the `SingleCells` and `MetaCells`
#' classes. den-SNE is t-SNE plus a density-preserving term: a tight population
#' stays tight in the embedding and a diffuse one stays diffuse. Plain t-SNE
#' inflates dense clusters and shrinks sparse ones, so relative cluster sizes
#' on the plot mean nothing; with den-SNE they do. Setting `lambda = 0` in
#' [manifoldsR::params_densne()] gives you back plain t-SNE.
#'
#' Neighbour handling, `approx_type` and the `k` versus `perplexity` caveats
#' are the same as in [bixverse::tsne_sc()].
#'
#' @param object `SingleCells`, `MetaCells` class.
#' @param use_knn Boolean. Use the kNN graph found in the object. Defaults to
#' `FALSE`. If not available, will default to the embedding.
#' @param embd_to_use String. The embedding to use for den-SNE. Must be
#' available in the object.
#' @param slot_name String. The name of this embedding within the object.
#' Defaults to `"densne"`.
#' @param no_embd_to_use Optional integer. Number of embedding dimensions to
#' use. If `NULL` all will be used.
#' @param modality String. On which modality to run den-SNE. One of
#' `c("rna", "adt", "wnn")`. The two latter options are only available for
#' multi-modal versions with the added data.
#' @param n_dim Integer. Number of den-SNE dimensions. Currently only `2L` is
#' supported. Defaults to `2L`.
#' @param perplexity Numeric. Perplexity parameter. Typical values between 5
#' and 50. Defaults to `10.0`.
#' @param approx_type String. Approximation method. One of `"bh"` (Barnes-Hut),
#' `"fft"` or `"fft_3k"`. Defaults to `"bh"`.
#' @param knn_method String. Approximate nearest neighbour algorithm. One of
#' `"hnsw"`, `"balltree"`, `"annoy"`, `"nndescent"`, or `"exhaustive"`.
#' @param nn_params Named list. See [manifoldsR::params_nn()].
#' @param tsne_params Named list. See [manifoldsR::params_tsne()].
#' @param dens_params Named list. The density knobs, see
#' [manifoldsR::params_densne()].
#' @param seed Integer. For reproducibility.
#' @param .verbose Boolean. Controls verbosity.
#'
#' @returns The object with a `"densne"` embedding added.
#'
#' @export
#'
#' @references Narayan, Berger & Cho, Nat. Biotechnol., 2021
#'
#' @examples
#' # Barnes-Hut den-SNE on the PCA factors
#' sc <- demo_single_cells()
#' sc <- densne_sc(sc, .verbose = FALSE)
#' dim(get_embedding(sc, "densne"))
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
densne_sc <- S7::new_generic(
  name = "densne_sc",
  dispatch_args = "object",
  fun = function(
    object,
    use_knn = FALSE,
    embd_to_use = "pca",
    slot_name = "densne",
    no_embd_to_use = NULL,
    modality = c("rna", "adt", "wnn"),
    n_dim = 2L,
    perplexity = 10.0,
    approx_type = c("bh", "fft", "fft_3k"),
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
    dens_params = manifoldsR::params_densne(),
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

S7::method(densne_sc, ScOrMc) <- function(
  object,
  use_knn = FALSE,
  embd_to_use = "pca",
  slot_name = "densne",
  no_embd_to_use = NULL,
  modality = c("rna", "adt", "wnn"),
  n_dim = 2L,
  perplexity = 10.0,
  approx_type = c("bh", "fft", "fft_3k"),
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
  dens_params = manifoldsR::params_densne(),
  seed = 42L,
  .verbose = TRUE
) {
  modality <- match.arg(modality)
  approx_type <- match.arg(approx_type)
  knn_method <- match.arg(knn_method)

  # checks
  checkmate::qassert(slot_name, "S1")
  checkmate::qassert(n_dim, "I1[2,2]")
  checkmate::qassert(perplexity, "N1[1,)")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, "B1")

  inputs <- .manifold_inputs(
    object = object,
    use_knn = use_knn,
    embd_to_use = embd_to_use,
    no_embd_to_use = no_embd_to_use,
    modality = modality
  )

  if (.verbose) {
    message("Running den-SNE.")
  }

  densne_embd <- manifoldsR::densne(
    data = inputs$embd,
    knn = inputs$knn,
    n_dim = n_dim,
    perplexity = perplexity,
    approx_type = approx_type,
    knn_method = knn_method,
    nn_params = nn_params,
    tsne_params = tsne_params,
    dens_params = dens_params,
    seed = seed,
    .verbose = .verbose
  )

  .store_manifold(
    object = object,
    embd = densne_embd,
    prefix = "densne",
    slot_name = slot_name,
    modality = modality,
    from = .manifold_from(
      embd_to_use = embd_to_use,
      cache_modality = inputs$cache_modality,
      modality = modality,
      has_knn = !is.null(inputs$knn)
    )
  )
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
#' @param modality String. On which modality to run PHATE. One of
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
#' @returns The object with a `"phate"` embedding added.
#'
#' @export
#'
#' @examples
#' # PHATE embedding off the cached kNN graph
#' sc <- demo_single_cells()
#' sc <- phate_sc(sc, .verbose = FALSE)
#' dim(get_embedding(sc, "phate"))
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
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
  knn_method <- match.arg(knn_method)

  # checks
  checkmate::qassert(slot_name, "S1")
  checkmate::qassert(n_dim, "I1[2,2]")
  checkmate::qassert(k, "I1[1,)")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, "B1")

  inputs <- .manifold_inputs(
    object = object,
    use_knn = use_knn,
    embd_to_use = embd_to_use,
    no_embd_to_use = no_embd_to_use,
    modality = modality
  )

  if (.verbose) {
    message("Running PHATE.")
  }

  phate_embd <- manifoldsR::phate(
    data = inputs$embd,
    knn = inputs$knn,
    n_dim = n_dim,
    k = k,
    knn_method = knn_method,
    nn_params = nn_params,
    phate_params = phate_params,
    seed = seed,
    .verbose = .verbose
  )

  .store_manifold(
    object = object,
    embd = phate_embd,
    prefix = "phate",
    slot_name = slot_name,
    modality = modality,
    from = .manifold_from(
      embd_to_use = embd_to_use,
      cache_modality = inputs$cache_modality,
      modality = modality,
      has_knn = !is.null(inputs$knn)
    )
  )
}

### forceatlas2 ----------------------------------------------------------------

#' Run ForceAtlas2 on a SingleCells/MetaCells object
#'
#' @description
#' Wrapper around [manifoldsR::forceatlas2()] and
#' [manifoldsR::forceatlas2_from_graph()] for the `SingleCells` and `MetaCells`
#' classes. ForceAtlas2 is a force-directed graph layout: every edge pulls,
#' every pair of cells pushes, and gravity keeps disconnected components from
#' drifting off. It is what scanpy's `draw_graph` does.
#'
#' `graph` picks what gets laid out:
#' \itemize{
#'   \item `"knn"` - The kNN graph (cached if `use_knn = TRUE`, otherwise built
#'   from the chosen embedding), turned into the UMAP fuzzy union graph first.
#'   This matches scanpy.
#'   \item `"snn"` - The sNN graph from [bixverse::find_neighbours_sc()], i.e.
#'   the same graph the Leiden/Louvain clustering runs on. `use_knn`,
#'   `embd_to_use`, `no_embd_to_use`, `k`, `knn_method` and `nn_params` do not
#'   apply here, and neither do the graph and initialisation knobs in
#'   `fa2_params`. Pass `init_embd` to start from an existing 2D embedding,
#'   otherwise the layout starts from random positions.
#' }
#'
#' @param object `SingleCells`, `MetaCells` class.
#' @param graph String. Which graph to lay out. One of `c("knn", "snn")`.
#' Defaults to `"knn"`.
#' @param use_knn Boolean. Use the kNN graph found in the object. Defaults to
#' `TRUE`. If not available, will default to the embedding. `"knn"` only.
#' @param embd_to_use String. The embedding to build the kNN graph from. Must
#' be available in the object. `"knn"` only.
#' @param slot_name String. The name of this embedding within the object.
#' Defaults to `"fa2"`.
#' @param no_embd_to_use Optional integer. Number of embedding dimensions to
#' use. If `NULL` all will be used. `"knn"` only.
#' @param init_embd Optional string. Name of a stored 2D embedding, e.g.
#' `"umap"`, to initialise the layout with. `"snn"` only. Defaults to `NULL`.
#' @param modality String. On which modality to run ForceAtlas2. One of
#' `c("rna", "adt", "wnn")`. The two latter options are only available for
#' multi-modal versions with the added data.
#' @param k Integer. Number of nearest neighbours. Defaults to `15L`. `"knn"`
#' only.
#' @param knn_method String. Approximate nearest neighbour algorithm. One of
#' `"hnsw"`, `"balltree"`, `"annoy"`, `"nndescent"`, or `"exhaustive"`.
#' `"knn"` only.
#' @param nn_params Named list. See [manifoldsR::params_nn()]. `"knn"` only.
#' @param fa2_params Named list. See [manifoldsR::params_fa2()].
#' @param seed Integer. For reproducibility.
#' @param .verbose Boolean. Controls verbosity.
#'
#' @returns The object with a `"fa2"` embedding added.
#'
#' @export
#'
#' @references Jacomy, et al., PLoS ONE, 2014
#'
#' @examples
#' # ForceAtlas2 on the cached kNN graph, then on the sNN graph
#' sc <- demo_single_cells()
#' sc <- forceatlas2_sc(sc, .verbose = FALSE)
#' sc <- forceatlas2_sc(
#'   sc,
#'   graph = "snn",
#'   slot_name = "fa2_snn",
#'   .verbose = FALSE
#' )
#' dim(get_embedding(sc, "fa2_snn"))
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
forceatlas2_sc <- S7::new_generic(
  name = "forceatlas2_sc",
  dispatch_args = "object",
  fun = function(
    object,
    graph = c("knn", "snn"),
    use_knn = TRUE,
    embd_to_use = "pca",
    slot_name = "fa2",
    no_embd_to_use = NULL,
    init_embd = NULL,
    modality = c("rna", "adt", "wnn"),
    k = 15L,
    knn_method = c(
      "kmknn",
      "hnsw",
      "balltree",
      "annoy",
      "nndescent",
      "exhaustive"
    ),
    nn_params = manifoldsR::params_nn(),
    fa2_params = manifoldsR::params_fa2(),
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

S7::method(forceatlas2_sc, ScOrMc) <- function(
  object,
  graph = c("knn", "snn"),
  use_knn = TRUE,
  embd_to_use = "pca",
  slot_name = "fa2",
  no_embd_to_use = NULL,
  init_embd = NULL,
  modality = c("rna", "adt", "wnn"),
  k = 15L,
  knn_method = c(
    "kmknn",
    "hnsw",
    "balltree",
    "annoy",
    "nndescent",
    "exhaustive"
  ),
  nn_params = manifoldsR::params_nn(),
  fa2_params = manifoldsR::params_fa2(),
  seed = 42L,
  .verbose = TRUE
) {
  graph <- match.arg(graph)
  modality <- match.arg(modality)
  knn_method <- match.arg(knn_method)

  # checks
  checkmate::qassert(slot_name, "S1")
  checkmate::qassert(init_embd, c("0", "S1"))
  checkmate::qassert(k, "I1[2,)")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, "B1")

  if (graph == "knn") {
    if (!is.null(init_embd)) {
      warning("`init_embd` only applies to graph = 'snn'. Ignoring it.")
    }

    inputs <- .manifold_inputs(
      object = object,
      use_knn = use_knn,
      embd_to_use = embd_to_use,
      no_embd_to_use = no_embd_to_use,
      modality = modality
    )

    if (.verbose) {
      message("Running ForceAtlas2 on the kNN graph.")
    }

    fa2_embd <- manifoldsR::forceatlas2(
      data = inputs$embd,
      knn = inputs$knn,
      k = k,
      knn_method = knn_method,
      nn_params = nn_params,
      fa2_params = fa2_params,
      seed = seed,
      .verbose = .verbose
    )

    from <- .manifold_from(
      embd_to_use = embd_to_use,
      cache_modality = inputs$cache_modality,
      modality = modality,
      has_knn = !is.null(inputs$knn)
    )
  } else {
    checkmate::assertTRUE(
      S7::S7_inherits(object, SingleCells) ||
        S7::S7_inherits(object, MetaCells) ||
        S7::S7_inherits(object, SingleCellsSubset)
    )
    if (modality != "rna" && !S7::S7_inherits(object, SingleCellsMultiModal)) {
      stop(sprintf(
        "modality = '%s' is only supported for SingleCellsMultiModal.",
        modality
      ))
    }

    # the graph and the init both live under the write modality, wnn included
    assert_sc_state(object, artefacts = "snn", modality = modality)
    snn_graph <- get_snn_graph(object, modality = modality)
    if (is.null(snn_graph)) {
      stop("No sNN graph found. Run find_neighbours_sc() first.")
    }

    init <- NULL
    if (!is.null(init_embd)) {
      assert_sc_state(object, artefacts = init_embd, modality = modality)
      init <- get_embedding(
        x = object,
        embd_name = init_embd,
        modality = modality
      )
    }

    if (.verbose) {
      message("Running ForceAtlas2 on the sNN graph.")
    }

    fa2_embd <- manifoldsR::forceatlas2_from_graph(
      graph = snn_graph,
      init = init,
      fa2_params = fa2_params,
      seed = seed,
      .verbose = .verbose
    )

    from <- c(
      sprintf("%s:snn", modality),
      if (!is.null(init_embd)) sprintf("%s:%s", modality, init_embd) else NULL
    )
  }

  .store_manifold(
    object = object,
    embd = fa2_embd,
    prefix = "fa2",
    slot_name = slot_name,
    modality = modality,
    from = from
  )
}

## feature extraction ----------------------------------------------------------

### helpers --------------------------------------------------------------------

#' Match requested features against available feature names
#'
#' @param features Character vector. Requested feature ids.
#' @param available Character vector. Feature names present in the data.
#'
#' @returns The subset of `features` that was matched, in input order.
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

#' Resolve feature names to gene indices, keeping the two in step
#'
#' @description
#' [get_gene_indices()] warns and silently drops anything it cannot match, so
#' the caller is left holding a `features` vector longer than the indices it
#' got back. Anything that labels a Rust result by position has to use the
#' surviving names rather than the requested ones, otherwise the labels slide
#' off the values. This is the in-memory `.match_features()` for the streaming
#' paths.
#'
#' @param object `SingleCells` or `SingleCellsSubset` class.
#' @param features Character vector. The requested gene ids.
#'
#' @returns A list with `features`, the surviving gene ids in request order,
#' and `indices`, their 0-indexed positions for Rust.
#'
#' @keywords internal
.resolve_gene_features <- function(object, features) {
  gene_idx <- get_gene_indices(
    x = object,
    gene_ids = features,
    rust_index = TRUE
  )

  list(
    features = unname(get_gene_names_from_idx(
      x = object,
      gene_idx = gene_idx,
      rust_based = TRUE
    )),
    indices = gene_idx
  )
}

#' Z-score and optionally clip a numeric vector
#'
#' @param x Numeric vector.
#' @param clip Optional numeric. Clip the z-scores to `[-clip, clip]`.
#'
#' @returns The scaled (and optionally clipped) vector.
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
#' @returns A data.table with a `cell_id` column and one column per feature.
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
#' @returns A long data.table with columns `gene`, `group`, `mean_exp`,
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
#' @returns The data.table with an added `scaled_exp` column and ordered
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

#' @method extract_dot_plot_data ScOrScSubset
S7::method(extract_dot_plot_data, ScOrScSubset) <- function(
  object,
  features,
  grouping_variable,
  scale_exp = TRUE,
  modality = c("rna", "adt")
) {
  modality <- match.arg(modality)
  checkmate::assertTRUE(
    S7::S7_inherits(object, SingleCells) ||
      S7::S7_inherits(object, SingleCellsSubset)
  )
  checkmate::qassert(features, "S+")
  checkmate::qassert(grouping_variable, "S1")
  checkmate::qassert(scale_exp, "B1")

  if (modality != "rna") {
    stop(paste(
      "SingleCells only supports modality = 'rna'.",
      "Use SingleCellsMultiModal for ADT."
    ))
  }

  resolved <- .resolve_gene_features(object, features)
  features <- resolved$features

  cell_idx <- get_cells_to_keep(object)
  grouping <- as.factor(
    unlist(object[[grouping_variable]], use.names = FALSE)
  )

  # NA group codes turn into an out-of-bounds usize Rust-side, so drop them
  keep <- !is.na(grouping)
  if (!any(keep)) {
    stop(sprintf("Grouping variable `%s` is entirely NA.", grouping_variable))
  }
  cell_idx <- cell_idx[keep]
  grouping <- droplevels(grouping[keep])

  gene_res <- rs_extract_grouped_gene_stats(
    f_path = get_rust_count_gene_f_path(object),
    cell_indices = cell_idx,
    gene_indices = resolved$indices,
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
    # lookup by concrete class: S7 expands unions at registration, but
    # S7::method() refuses to retrieve by one
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

#' @method extract_gene_expression ScOrScSubset
S7::method(extract_gene_expression, ScOrScSubset) <- function(
  object,
  features,
  obs_cols = NULL,
  scale = FALSE,
  clip = NULL,
  modality = c("rna", "adt"),
  layer = c("norm", "magic")
) {
  modality <- match.arg(modality)
  layer <- match.arg(layer)
  checkmate::assertTRUE(
    S7::S7_inherits(object, SingleCells) ||
      S7::S7_inherits(object, SingleCellsSubset)
  )
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

  dt <- if (identical(layer, "magic")) {
    # the imputed layer is already dense and in memory, so it takes the same
    # route as the ADT counts do
    mat <- .magic_matrix(object)
    .extract_expr_in_memory(
      mat,
      .match_features(features, colnames(mat)),
      scale,
      clip
    )
  } else {
    resolved <- .resolve_gene_features(object, features)

    counts <- rs_extract_several_genes_plots(
      f_path = get_rust_count_gene_f_path(object),
      cell_indices = get_cells_to_keep(object),
      gene_indices = resolved$indices,
      scale = scale,
      clip = clip
    )

    out <- data.table::data.table(
      cell_id = get_cell_names(object, filtered = TRUE)
    )
    for (i in seq_along(resolved$features)) {
      data.table::set(out, j = resolved$features[i], value = counts[[i]])
    }
    out
  }

  if (!is.null(obs_cols)) {
    obs_dt <- object[[obs_cols]]
    for (col in names(obs_dt)) {
      data.table::set(dt, j = col, value = obs_dt[[col]])
    }
  }

  dt
}

#' Fetch the imputed matrix, with a pointer when it is missing
#'
#' @param object `SingleCells` or `SingleCellsSubset` class.
#'
#' @returns The cells x genes matrix held by the `ScMagic` layer.
#'
#' @keywords internal
.magic_matrix <- function(object) {
  magic <- get_magic(object)

  if (is.null(magic)) {
    stop(paste(
      "No imputed layer found. Run run_magic_sc() first,",
      "or use layer = 'norm'."
    ))
  }

  magic[["data"]]
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
  modality = c("rna", "adt"),
  layer = c("norm", "magic")
) {
  modality <- match.arg(modality)
  layer <- match.arg(layer)

  if (modality == "rna") {
    # lookup by concrete class: S7 expands unions at registration, but
    # S7::method() refuses to retrieve by one
    rna_method <- S7::method(extract_gene_expression, SingleCells)
    return(rna_method(
      object = object,
      features = features,
      obs_cols = obs_cols,
      scale = scale,
      clip = clip,
      modality = "rna",
      layer = layer
    ))
  }

  # ADT path
  if (identical(layer, "magic")) {
    stop(paste(
      "The imputed layer holds RNA counts.",
      "Use modality = 'rna' or layer = 'norm'."
    ))
  }

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
  modality = c("rna", "adt"),
  layer = c("norm", "magic")
) {
  modality <- match.arg(modality)
  layer <- match.arg(layer)
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

  if (identical(layer, "magic")) {
    stop(paste(
      "MAGIC is not available for MetaCells, which are already aggregated.",
      "Use layer = 'norm'."
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
