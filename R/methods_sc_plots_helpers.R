# single cell plot helpers -----------------------------------------------------

## 2d embeddings ---------------------------------------------------------------

### umap -----------------------------------------------------------------------

#' Run UMAP on a SingleCells object
#'
#' @description
#' Wrapper around [manifoldsR::umap()] for the `SingleCells` class.
#'
#' @param object `SingleCells` class.
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
#' `"balltree"`, `"hnsw"`, `"annoy"`, `"nndescent"`, or `"exhaustive"`.
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
    embd_to_use = "pca",
    no_embd_to_use = NULL,
    n_dim = 2L,
    k = 15L,
    min_dist = 0.5,
    spread = 1.0,
    knn_method = c("balltree", "hnsw", "annoy", "nndescent", "exhaustive"),
    nn_params = manifoldsR::params_nn(),
    umap_params = manifoldsR::params_umap(),
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method umap_sc SingleCells
#'
#' @export
S7::method(umap_sc, SingleCells) <- function(
  object,
  embd_to_use = "pca",
  no_embd_to_use = NULL,
  n_dim = 2L,
  k = 15L,
  min_dist = 0.5,
  spread = 1.0,
  knn_method = c("balltree", "hnsw", "annoy", "nndescent", "exhaustive"),
  nn_params = manifoldsR::params_nn(),
  umap_params = manifoldsR::params_umap(),
  seed = 42L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::qassert(embd_to_use, "S1")
  checkmate::qassert(no_embd_to_use, c("I1", "0"))
  checkmate::qassert(n_dim, "I1[1,)")
  checkmate::qassert(k, "I1[2,)")
  checkmate::qassert(min_dist, "N1[0,)")
  checkmate::qassert(spread, "N1[0,)")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, "B1")
  knn_method <- match.arg(knn_method)

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
#' Wrapper around [manifoldsR::tsne()] for the `SingleCells` class.
#'
#' @param object `SingleCells` class.
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
#' `"balltree"`, `"hnsw"`, `"annoy"`, `"nndescent"`, or `"exhaustive"`.
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
    embd_to_use = "pca",
    no_embd_to_use = NULL,
    n_dim = 2L,
    perplexity = 30.0,
    approx_type = c("bh", "fft"),
    knn_method = c("balltree", "hnsw", "annoy", "nndescent", "exhaustive"),
    nn_params = manifoldsR::params_nn(),
    tsne_params = manifoldsR::params_tsne(),
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method tsne_sc SingleCells
#'
#' @export
S7::method(tsne_sc, SingleCells) <- function(
  object,
  embd_to_use = "pca",
  no_embd_to_use = NULL,
  n_dim = 2L,
  perplexity = 30.0,
  approx_type = c("bh", "fft"),
  knn_method = c("balltree", "hnsw", "annoy", "nndescent", "exhaustive"),
  nn_params = manifoldsR::params_nn(),
  tsne_params = manifoldsR::params_tsne(),
  seed = 42L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::qassert(embd_to_use, "S1")
  checkmate::qassert(no_embd_to_use, c("I1", "0"))
  checkmate::qassert(n_dim, "I1[2,2]")
  checkmate::qassert(perplexity, "N1[1,)")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, "B1")
  approx_type <- match.arg(approx_type)
  knn_method <- match.arg(knn_method)

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
#' Wrapper around [manifoldsR::phate()] for the `SingleCells` class.
#'
#' @param object `SingleCells` class.
#' @param embd_to_use String. The embedding to use for PHATE. Must be available
#' in the object.
#' @param no_embd_to_use Optional integer. Number of embedding dimensions to
#' use. If `NULL` all will be used.
#' @param n_dim Integer. Number of PHATE dimensions. Currently only `2L` is
#' supported. Defaults to `2L`.
#' @param k Integer. Number of nearest neighbours for graph construction.
#' Defaults to `5L`.
#' @param knn_method String. Approximate nearest neighbour algorithm. One of
#' `"balltree"`, `"hnsw"`, `"annoy"`, `"nndescent"`, or `"exhaustive"`.
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
    embd_to_use = "pca",
    no_embd_to_use = NULL,
    n_dim = 2L,
    k = 5L,
    knn_method = c("balltree", "hnsw", "annoy", "nndescent", "exhaustive"),
    nn_params = manifoldsR::params_nn(),
    phate_params = manifoldsR::params_phate(),
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method phate_sc SingleCells
#'
#' @export
S7::method(phate_sc, SingleCells) <- function(
  object,
  embd_to_use = "pca",
  no_embd_to_use = NULL,
  n_dim = 2L,
  k = 5L,
  knn_method = c("balltree", "hnsw", "annoy", "nndescent", "exhaustive"),
  nn_params = manifoldsR::params_nn(),
  phate_params = manifoldsR::params_phate(),
  seed = 42L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::qassert(embd_to_use, "S1")
  checkmate::qassert(no_embd_to_use, c("I1", "0"))
  checkmate::qassert(n_dim, "I1[2,2]")
  checkmate::qassert(k, "I1[1,)")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, "B1")
  knn_method <- match.arg(knn_method)

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
