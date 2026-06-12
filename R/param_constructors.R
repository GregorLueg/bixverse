# param constructors -----------------------------------------------------------

## defaults --------------------------------------------------------------------

#' Helper function to generate kNN defaults
#'
#' @description
#' This function generates various sensible default parameters for all of the
#' different approximate nearest neighbours that are available within this
#' package.
#'
#' @returns A list with default parameters for kNN searches. Following
#' parameters:
#' \itemize{
#'  \item k - Number of neighbours. Defaults to `15L`.
#'  \item knn_method - Which of method to use for the approximate nearest
#'  neighbour search. Defaults to `"kmknn"`. The implementations are:
#'  `c("kmknn", "hnsw", "annoy", "nndescent", "ivf", "exhaustive")`.
#'  \item ann_dist - Which distance metric to use for the approximate nearest
#'  neighbour search. Defaults to `"cosine"`. The implementations are
#'  `c("cosine", "euclidean")`.
#'  \item n_trees - Annoy param: number of trees to generate for Annoy. Defaults
#'  to `50L`.
#'  \item search_budget - Annoy param: optional search budget per tree for
#'  Annoy. If not provided, it will default to `n_tree * k * 20L`.
#'  \item diversify_prob - NNDescent param: diversification probability for the
#'  NNDescent index. This will diversify the index at the end and identify
#'  potentiall better edges. Defaults to `0.0`.
#'  \item delta - NNDescent param: early termination criterium for NNDescent.
#'  Defaults to `0.001`.
#'  \item ef_budget - NNDescent param: optional query budget parameter. Can
#'  accelerate querying, but at the cost of Recall.
#'  \item m - HNSW param: number of connections between layers for HNSW.
#'  Defaults to `16L`.
#'  \item ef_construction - HNSW param: size of dynamic candidate list during
#'  construction. Defaults to `200L`.
#'  \item ef_search - HNSW param: size of candidate list (higher = better
#'  recall, slower). Defaults to `100L`.
#'  \item n_list - IVF param: number of clusters/centroids to generate. Defaults
#'  to `NULL` (sqrt(n) n_lists will be generated in this case).
#'  \item n_probe - IVF param: number of clusters/centroids to query Defaults
#'  to `NULL` (sqrt(n_lists) clusters will be queried in this case).
#' }
#'
#' @export
params_knn_defaults <- function() {
  list(
    # General parameters
    k = 15L,
    knn_method = "kmknn",
    ann_dist = "euclidean",
    # Annoy
    n_trees = 50L,
    search_budget = NULL,
    # NNDescent
    delta = 0.001,
    diversify_prob = 0.0,
    ef_budget = NULL,
    # HNSW
    m = 16L,
    ef_construction = 200L,
    ef_search = 100L,
    # IVF
    n_list = NULL,
    n_probe = NULL
  )
}

#' Helper function to generate HVG defaults
#'
#' @returns A list with default parameters for kNN searches. Following
#' parameters:
#' \itemize{
#'  \item min_gene_var_pctl - Which percentile of the highly variable genes
#'  to include. Defaults to `0.7`.
#'  \item hvg_method - Which method to use to identify HVG. Defaults to `"vst"`.
#'  \item loess_span - In case of `"vst"` the span of the loess function.
#'  \item clip_max - The maximum clipping value (optional).
#'  \item n_bins - The number of bins to use for the `"meanvarbin"` HVG
#'  detection.
#'  \item binning_strategy - Which binning strategy to use for `"meanvarbin"`.
#' }
#'
#' @export
params_hvg_defaults <- function() {
  list(
    min_gene_var_pctl = 0.7,
    hvg_method = "vst",
    loess_span = 0.3,
    clip_max = NULL,
    n_bins = 20L,
    binning_strategy = "equal_width"
  )
}

#' Helper function to generate normalisation defaults for doublet detection.
#'
#' @return A list with the following parameters for normalisation specifically
#' designed for doublet detection methods:
#' \itemize{
#'  \item log_transform - Boolean. Shall the counts be log-normalised.
#'  Defaults to `TRUE`.
#'  \item mean_center - Boolean. Shall mean centring be applied. Defaults
#'  to `FALSE`.
#'  \item normalise_variance - Boolean. Shall the variance be normalised.
#'  Defaults to `FALSE`.
#'  \item target_size - Target library size. Defaults to `1e6`
#' }
#'
#' @export
params_norm_doublets_defaults <- function() {
  list(
    log_transform = TRUE,
    mean_center = FALSE,
    normalise_variance = FALSE,
    target_size = 1e6
  )
}

#' Helper function to generate default parameters for PCA
#'
#' @return A list with the following parameters for PCA.
#' \itemize{
#'  \item no_pcs - Integer. Number of PCs to consider. Defaults to `30L`.
#'  \item random_svd - Boolean. Shall randomised SVD be used. Defaults to
#'  `TRUE`.
#'  \item sparse - Boolean. Shall sparse solvers be used that do not do
#'  scaling. If set to yes, in the case of `random_svd = FALSE`, Lanczos
#'  iterations are used to solve the sparse SVD. With `random_svd = TRUE`, the
#'  sparse initial matrix is multiplied with the random matrix, yielding a
#'  much smaller dense matrix that does not increase the memory pressure
#'  massively.
#' }
params_pca_defaults <- function() {
  list(
    no_pcs = 30L,
    random_svd = TRUE,
    sparse = FALSE
  )
}

#' Helper function to generate default parameters for the fast clustering for
#' the doublet detection methods
#'
#' @returns A list with the following parameters for fast clustering:
#' \itemize{
#'  \item km_type - The type of k-means clustering. Defaults to `"minibatch"`
#'  \item n_centroids - The number of centroids to use. Default to `NULL` and
#'  the function will use `sqrt(N_cells) * 4` for the number of n_centroids.
#'  \item kmeans_iters - Number of maximum k-means iterations. Defaults to
#'  `100L`
#'  \item batch_size - Max batch size will be set to `4098L`, but pending data
#'  set set to `N_cells / 2`.
#' }
params_fast_cluster_default <- function() {
  list(
    km_type = "minibatch",
    n_centroids = NULL,
    kmeans_iters = 100L,
    batch_size = 4098L
  )
}

#' K-mean parameter defaults.
#'
#' @description
#' Helper function to generate defaults for the k-mean clustering were more
#' control is needed.
#'
#' @returns A list with the following parameters
#' \itemize{
#'  \item k_means_iter - Integer. The number of iterations to use for the
#'  clustering.
#'  \item k_means_init - String. The initialisation. Options are `"random"` and
#'  `"parallel"`. Defaults to `"parallel"`.
#'  \item gemm - Optional boolean. Controls which CPU implementation is used
#'  by the method. GEMM is faster with large dimensionality.
#'  \item hamerly - Optional boolean. Shall a faster exact method be used
#'  leveraging the triangle inequality. Faster on large data sets with large
#'  numbers of centroids.
#' }
params_kmeans_defaults <- function() {
  list(
    k_means_iter = 30L,
    k_means_init = "parallel",
    gemm = FALSE,
    hamerly = TRUE
  )
}

# constructors -----------------------------------------------------------------

## doublet detections ----------------------------------------------------------

### scrublet -------------------------------------------------------------------

#' Wrapper function for Scrublet doublet detection parameters
#'
#' @description Constructor for the various Scrublet parameters. In this case,
#' the default for the kNN graph generation was set to `"hnsw"` as this
#' algorithm showed the best performance in different empirical benchmarks.
#'
#' @param sim_doublet_ratio Numeric. Number of doublets to simulate relative to
#' the number of observed cells. For example, 2.0 simulates twice as many
#' doublets as there are cells. Defaults to `1.5`.
#' @param expected_doublet_rate Numeric. Expected doublet rate for the
#' experiment, typically 0.05-0.10 depending on cell loading. Must be between
#' 0 and 1. Defaults to `0.1`.
#' @param stdev_doublet_rate Numeric. Uncertainty in the expected doublet rate.
#' Defaults to `0.02`.
#' @param n_bins_histogram Integer. Number of bins for histogram-based automatic
#' threshold detection. Typically 50-100. Defaults to `100L`.
#' @param manual_threshold Optional numeric. Manual doublet score threshold. If
#' `NULL` (default), threshold is automatically detected from simulated doublet
#' score distribution.
#' @param normalisation List. Optional overrides for normalisation parameters.
#' See [bixverse::params_norm_doublets_defaults()] for available
#' parameters: `log_transform`, `mean_center`, `normalise_variance`,
#' `target_size`.
#' @param hvg List. Optional overrides for highly variable gene selection
#' parameters. See [bixverse::params_hvg_defaults()] for available parameters:
#' `min_gene_var_pctl`, `hvg_method`, `loess_span`, `clip_max`.
#' @param pca List. Optional overrides for PCA parameters. See
#' [bixverse::params_pca_defaults()] for available parameters: `no_pcs`,
#' `random_svd`, `sparse` and `skip_first_pc`.
#' @param knn List. Optional overrides for kNN parameters. See
#' [bixverse::params_knn_defaults()] for available parameters: `k`,
#' `knn_method`, `ann_dist`, `search_budget`, `n_trees`, `delta`,
#' `diversify_prob`, `ef_budget`, `m`, `ef_construction`, `ef_search`, `n_list`
#' and `n_probe`.
#' Note: this function defaults to `k = 0L` (automatic neighbour detection).
#'
#' @returns A named list with all Scrublet parameters, combining defaults with
#' any user-specified overrides.
#'
#' @export
params_scrublet <- function(
  sim_doublet_ratio = 1.5,
  expected_doublet_rate = 0.1,
  stdev_doublet_rate = 0.02,
  n_bins_histogram = 100L,
  manual_threshold = NULL,
  normalisation = list(),
  hvg = list(),
  pca = list(),
  knn = list(k = 0L)
) {
  # doublet simulation checks
  checkmate::qassert(sim_doublet_ratio, "N1(0,)")
  checkmate::qassert(expected_doublet_rate, "N1[0,1]")
  checkmate::qassert(stdev_doublet_rate, "N1[0,1]")
  checkmate::qassert(n_bins_histogram, "I1[10,)")
  if (!is.null(manual_threshold)) {
    checkmate::qassert(manual_threshold, "N1[0,)")
  }

  # generate final parameters
  params <- list(
    normalisation = modifyList(
      params_norm_doublets_defaults(),
      normalisation,
      keep.null = TRUE
    ),
    hvg = modifyList(params_hvg_defaults(), hvg, keep.null = TRUE),
    pca = modifyList(params_pca_defaults(), pca, keep.null = TRUE),
    knn = modifyList(params_knn_defaults(), knn, keep.null = TRUE),
    sim_doublet_ratio = sim_doublet_ratio,
    expected_doublet_rate = expected_doublet_rate,
    stdev_doublet_rate = stdev_doublet_rate,
    n_bins_histogram = n_bins_histogram,
    manual_threshold = manual_threshold
  )

  params <- purrr::list_flatten(params, name_spec = "{inner}")

  params
}

### boost ----------------------------------------------------------------------

#' Wrapper function for Boost parameters
#'
#' @param boost_rate Numeric. Boosting rate for the algorithm. Must be between
#' 0 and 1. Defaults to `0.25`.
#' @param replace Boolean. Whether to use replacement during boosting. Defaults
#' to `FALSE`.
#' @param resolution Numeric. Resolution parameter for graph-based clustering.
#' Higher values lead to more clusters. Defaults to `1.0`.
#' @param n_iters Integer. Number of iterations to run the algorithm. Defaults
#' to `20L`.
#' @param p_thresh Numeric. P-value threshold for significance testing. Defaults
#' to `1e-7`.
#' @param voter_thresh Numeric. Voter threshold across iterations. Proportion of
#' iterations a cell must be assigned to a cluster to be considered a member.
#' Must be between 0 and 1. Defaults to `0.9`.
#' @param fast_cluster Boolean. Shall fast Louvain clustering be applied, i.e.,
#' k-means clustering and use the centroids for kNN graph generation and
#' Louvain clustering with then backpropagating the membership based on centroid
#' proximity.
#' @param normalisation List. Optional overrides for normalisation parameters.
#' See [bixverse::params_norm_doublets_defaults()] for available
#' parameters: `log_transform`, `mean_center`, `normalise_variance`,
#' `target_size`. Note: Boost uses different defaults (`log_transform = FALSE`,
#' `mean_center = TRUE`, `normalise_variance = TRUE`, `target_size = NULL`).
#' @param hvg List. Optional overrides for highly variable gene selection
#' parameters. See [bixverse::params_hvg_defaults()] for available parameters:
#' `min_gene_var_pctl`, `hvg_method`, `loess_span`, `clip_max`.
#' @param pca List. Optional overrides for PCA parameters. See
#' [bixverse::params_pca_defaults()] for available parameters: `no_pcs`,
#' `random_svd`.
#' @param knn List. Optional overrides for kNN parameters. See
#' [bixverse::params_knn_defaults()] for available parameters: `k`,
#' `knn_method`, `ann_dist`, `search_budget`, `n_trees`, `delta`,
#' `diversify_prob`, `ef_budget`, `m`, `ef_construction`, `ef_search`, `n_list`
#' and `n_probe`. Note: this function defaults to `k = 0L` (automatic neighbour
#' detection).
#' @param fast_cluster_params List. Optional overrides for the fast clustering
#' parameters. Only relevant if `fast_cluster = TRUE`. See
#' [params_fast_cluster_default()] for available parameters: `km_type`,
#' `n_centroids`, `kmeans_iters` and `batch_size`.
#'
#' @returns A named list with all Boost parameters, combining defaults with
#' any user-specified overrides.
#'
#' @returns A list with the Boost parameters.
#'
#' @export
params_boost <- function(
  boost_rate = 0.25,
  replace = FALSE,
  resolution = 1.0,
  n_iters = 20L,
  p_thresh = 1e-7,
  voter_thresh = 0.9,
  fast_cluster = FALSE,
  normalisation = list(),
  hvg = list(),
  pca = list(),
  knn = list(k = 0L),
  fast_cluster_params = list()
) {
  # checks
  checkmate::qassert(boost_rate, "N1[0,1]")
  checkmate::qassert(replace, "B1")
  checkmate::qassert(resolution, "N1(0,)")
  checkmate::qassert(n_iters, "I1[1,)")
  checkmate::qassert(p_thresh, "N1(0,)")
  checkmate::qassert(voter_thresh, "N1[0,1]")
  checkmate::qassert(fast_cluster, "B1")

  # generate final parameters
  params <- list(
    normalisation = modifyList(
      params_norm_doublets_defaults(),
      normalisation,
      keep.null = TRUE
    ),
    hvg = modifyList(params_hvg_defaults(), hvg, keep.null = TRUE),
    pca = modifyList(params_pca_defaults(), pca, keep.null = TRUE),
    knn = modifyList(params_knn_defaults(), knn, keep.null = TRUE),
    fast_cluster_params = modifyList(
      params_fast_cluster_default(),
      fast_cluster_params,
      keep.null = TRUE
    ),
    boost_rate = boost_rate,
    replace = replace,
    resolution = resolution,
    fast_cluster = fast_cluster,
    n_iters = n_iters,
    p_thresh = p_thresh,
    voter_thresh = voter_thresh
  )

  params <- purrr::list_flatten(params, name_spec = "{inner}")

  assertScBoost(params)

  params
}

### scdblfinder ----------------------------------------------------------------

#' Wrapper function for scDblFinder doublet detection parameters
#'
#' @description Constructor for the scDblFinder parameters. This method
#' combines cluster-aware doublet simulation with a gradient-boosted
#' classifier trained on engineered features.
#'
#' @param n_genes Integer. Number of top-expressed genes to use as features.
#' Defaults to `1352L`.
#' @param doublet_ratio Numeric. Ratio of simulated doublets to observed cells.
#' Defaults to `1.0`.
#' @param heterotypic_bias Numeric. Fraction of simulated pairs forced to come
#' from different clusters (0-1). Defaults to `1.0`.
#' @param cluster_resolution Numeric. Resolution for the initial Louvain
#' clustering. Defaults to `1.0`.
#' @param cluster_iters Integer. Number of Louvain iterations per clustering
#' step. Defaults to `10L`.
#' @param fast_cluster Boolean. Shall fast Louvain clustering be applied, i.e.,
#' k-means clustering and use the centroids for kNN graph generation and
#' Louvain clustering with then backpropagating the membership based on centroid
#' proximity.
#' @param n_iterations Integer. Number of refinement iterations. Typically 2-3.
#' Defaults to `3L`.
#' @param n_trees Integer. Maximum number of boosting rounds for the GBM
#' classifier. Defaults to `200L`.
#' @param max_depth Integer. Maximum tree depth. Shallow trees (3-5) work best.
#' Defaults to `4L`.
#' @param learning_rate Numeric. Shrinkage applied to each tree. Defaults to
#' `0.3`.
#' @param min_samples_leaf Integer. Minimum training samples per leaf. Defaults
#' to `20L`.
#' @param subsample_rate Numeric. Fraction of samples used per tree. Defaults to
#' `0.75`.
#' @param cv_folds Integer. Number of cross-validation folds for boosting round
#' selection. Defaults to `5L`.
#' @param cv_early_stop Integer. Early stopping patience per CV fold. Defaults
#' to `2L`.
#' @param se_fraction Numeric. Multiplier on the standard error for the SE rule
#' used in round selection. Defaults to `1.0`
#' @param include_pcs Integer. Number of leading principal components to include
#' as classifier features. Defaults to `19L`.
#' @param expected_doublet_rate Optional numeric. Expected doublet rate as a
#' percentage. If not provided, will be calculated internally.
#' @param manual_threshold Optional numeric. Manual score threshold. If `NULL`
#' (default), expected-rate thresholding is used.
#' @param normalisation List. Optional overrides for normalisation parameters.
#' See [bixverse::params_norm_doublets_defaults()].
#' @param pca List. Optional overrides for PCA parameters.
#' See [bixverse::params_pca_defaults()].
#' @param knn List. Optional overrides for kNN parameters.
#' See [bixverse::params_knn_defaults()]. NNDescent works better for the larger
#' k-values often used here.
#' @param fast_cluster_params List. Optional overrides for the fast clustering
#' parameters. Only relevant if `fast_cluster = TRUE`. See
#' [params_fast_cluster_default()] for available parameters: `km_type`,
#' `n_centroids`, `kmeans_iters` and `batch_size`.
#'
#' @returns A named list with all scDblFinder parameters.
#'
#' @export
params_scdblfinder <- function(
  n_genes = 1352L,
  doublet_ratio = 1.0,
  heterotypic_bias = 1.0,
  cluster_resolution = 1.0,
  cluster_iters = 10L,
  fast_cluster = FALSE,
  n_iterations = 3L,
  n_trees = 200L,
  max_depth = 4L,
  learning_rate = 0.3,
  min_samples_leaf = 20L,
  subsample_rate = 0.75,
  cv_folds = 5L,
  cv_early_stop = 2L,
  se_fraction = 1.0,
  include_pcs = 19L,
  expected_doublet_rate = NULL,
  cxds_genes = NULL,
  manual_threshold = NULL,
  normalisation = list(mean_center = TRUE),
  pca = list(),
  knn = list(k = 0L),
  fast_cluster_params = list()
) {
  # checks
  checkmate::qassert(n_genes, "I1[1,)")
  checkmate::qassert(doublet_ratio, "N1(0,)")
  checkmate::qassert(heterotypic_bias, "N1[0,1]")
  checkmate::qassert(cluster_resolution, "N1(0,)")
  checkmate::qassert(cluster_iters, "I1[1,)")
  checkmate::qassert(n_iterations, "I1[1,)")
  checkmate::qassert(fast_cluster, "B1")
  checkmate::qassert(n_trees, "I1[1,)")
  checkmate::qassert(max_depth, "I1[1,)")
  checkmate::qassert(learning_rate, "N1(0,)")
  checkmate::qassert(min_samples_leaf, "I1[1,)")
  checkmate::qassert(subsample_rate, "N1(0,1]")
  checkmate::qassert(cv_folds, "I1[2,)")
  checkmate::qassert(cv_early_stop, "I1[1,)")
  checkmate::qassert(se_fraction, "N1[0,)")
  checkmate::qassert(expected_doublet_rate, c("N1(0,1]", "0"))
  checkmate::qassert(cxds_genes, c("I1", "0"))
  checkmate::qassert(manual_threshold, c("N1[0,)", "0"))

  # generate params list
  params <- list(
    normalisation = modifyList(
      params_norm_doublets_defaults(),
      normalisation,
      keep.null = TRUE
    ),
    pca = modifyList(params_pca_defaults(), pca, keep.null = TRUE),
    knn = modifyList(params_knn_defaults(), knn, keep.null = TRUE),
    fast_cluster_params = modifyList(
      params_fast_cluster_default(),
      fast_cluster_params,
      keep.null = TRUE
    ),
    n_genes = n_genes,
    doublet_ratio = doublet_ratio,
    heterotypic_bias = heterotypic_bias,
    cluster_resolution = cluster_resolution,
    cluster_iters = cluster_iters,
    fast_cluster = fast_cluster,
    n_iterations = n_iterations,
    n_trees = n_trees,
    max_depth = max_depth,
    learning_rate = learning_rate,
    min_samples_leaf = min_samples_leaf,
    subsample_rate = subsample_rate,
    cv_folds = cv_folds,
    cv_early_stop = cv_early_stop,
    se_fraction = se_fraction,
    include_pcs = include_pcs,
    expected_doublet_rate = expected_doublet_rate,
    manual_threshold = manual_threshold,
    cxds_genes = cxds_genes
  )

  params <- purrr::list_flatten(params, name_spec = "{inner}")

  params
}

## neighbours ------------------------------------------------------------------

#' Wrapper function for parameters for neighbour identification in single cell
#'
#' @param full_snn Boolean. Shall the full shared nearest neighbour graph
#' be generated that generates edges between all cells instead of between
#' only neighbours.
#' @param pruning Numeric. Weights below this threshold will be set to 0 in
#' the generation of the sNN graph. Seurat uses for example `1/15` with
#' `k = 20`. As the default k is set to 15, we set it to `1/12`.
#' @param snn_similarity String. One of `c("rank", "jaccard")`. The Jaccard
#' similarity calculates the Jaccard between the neighbours, whereas the rank
#' method calculates edge weights based on the ranking of shared neighbours.
#' For the rank method, the weight is determined by finding the shared
#' neighbour with the lowest combined rank across both cells, where
#' lower-ranked (closer) shared neighbours result in higher edge weights
#' Both methods produce weights normalised to the range `[0, 1]`.
#' @param knn List. Optional overrides for kNN parameters. See
#' [bixverse::params_knn_defaults()] for available parameters: `k`,
#' `knn_method`, `ann_dist`, `search_budget`, `n_trees`, `delta`,
#' `diversify_prob`, `ef_budget`, `m`, `ef_construction`, `ef_search`, `n_list`
#' and `n_probe`.
#'
#' @returns A list with the neighbour parameters.
#'
#' @export
params_sc_neighbours <- function(
  full_snn = TRUE,
  pruning = 1 / 12,
  snn_similarity = c("jaccard", "rank"),
  knn = list(ann_dist = "cosine")
) {
  snn_similarity <- match.arg(snn_similarity)

  checkmate::qassert(full_snn, "B1")
  checkmate::qassert(pruning, "N1[0, 1]")
  checkmate::assertChoice(snn_similarity, c("rank", "jaccard"))

  knn_params <- modifyList(
    params_knn_defaults(),
    knn,
    keep.null = TRUE
  )

  c(
    list(
      full_snn = full_snn,
      pruning = pruning,
      snn_similarity = snn_similarity
    ),
    knn_params
  )
}

## fast clustering -------------------------------------------------------------

#' Fast single cell clustering parameters
#'
#' @param kmeans_iters Integer. Number of iterations for k-means clustering.
#' @param batch_size Integer. Batch size for mini batch k-means clustering.
#' @param drift_threshold Numeric. The drift for the mini batch k-means
#' clustering. If the centroid drift is below this, the mini batch k-means
#' terminates.
#' @param lr_alpha Numeric. Learning rate alpha parameter for mini batch
#' k-means.
#' @param louvain_iters Integer. Number of iterations for Louvain clustering.
#' @param full_snn Boolean. Shall the full shared nearest neighbour graph
#' be generated that generates edges between all cells instead of between
#' only neighbours.
#' @param pruning Optional numeric. Weights below this threshold will be set to
#' 0 in the generation of the sNN graph. If not provided, defaults to
#' `1 / ceil(k * 0.8)`.
#' @param snn_similarity String. One of `c("rank", "jaccard")`. The Jaccard
#' similarity calculates the Jaccard between the neighbours, whereas the rank
#' method calculates edge weights based on the ranking of shared neighbours.
#' For the rank method, the weight is determined by finding the shared
#' neighbour with the lowest combined rank across both cells, where
#' lower-ranked (closer) shared neighbours result in higher edge weights
#' Both methods produce weights normalised to the range `[0, 1]`.
#' @param knn List. Optional overrides for kNN parameters. See
#' [bixverse::params_knn_defaults()] for available parameters: `k`,
#' `knn_method`, `ann_dist`, `search_budget`, `n_trees`, `delta`,
#' `diversify_prob`, `ef_budget`, `m`, `ef_construction`, `ef_search`, `n_list`
#' and `n_probe`. Sets the default `k = 5L`.
#'
#' @returns A named list with the single cell fast clustering parameters.
#'
#' @export
params_sc_fast_cluster <- function(
  # kmeans
  kmeans_iters = 100L,
  batch_size = 4096L,
  drift_threshold = 1e-4,
  lr_alpha = 1.0,
  # snn
  full_snn = FALSE,
  pruning = NULL,
  snn_similarity = c("jaccard", "rank"),
  # louvain
  louvain_iters = 10L,
  # knn
  knn = list(k = 5L)
) {
  snn_similarity <- match.arg(snn_similarity)

  # checks
  checkmate::qassert(kmeans_iters, "I1")
  checkmate::qassert(batch_size, "I1")
  checkmate::qassert(drift_threshold, "N1")
  checkmate::qassert(lr_alpha, "N1")
  checkmate::qassert(louvain_iters, "I1")
  checkmate::qassert(full_snn, "B1")
  checkmate::qassert(pruning, c("N1", "0"))
  checkmate::assertChoice(snn_similarity, c("jaccard", "rank"))

  knn_params <- modifyList(
    params_knn_defaults(),
    knn,
    keep.null = TRUE
  )

  c(
    list(
      kmeans_iters = kmeans_iters,
      batch_size = batch_size,
      drift_threshold = drift_threshold,
      lr_alpha = lr_alpha,
      louvain_iters = louvain_iters,
      full_snn = full_snn,
      pruning = pruning,
      snn_similarity = snn_similarity
    ),
    knn_params
  )
}

## vision ----------------------------------------------------------------------

#' Wrapper function for parameters for VISION with auto-correlation
#'
#' @param n_perm Integer. Number of random gene sets to generate per cluster.
#' @param n_cluster Integer. Number of clusters for the random gene set
#' clustering generation.
#' @param knn List. Optional overrides for kNN parameters. See
#' [bixverse::params_knn_defaults()] for available parameters: `k`,
#' `knn_method`, `ann_dist`, `search_budget`, `n_trees`, `delta`,
#' `diversify_prob`, `ef_budget`, `m`, `ef_construction`, `ef_search`, `n_list`
#' and `n_probe`.
#'
#' @returns A list with the VISION parameters when you wish to use the
#' auto-correlation version.
#'
#' @export
params_sc_vision <- function(
  n_perm = 500L,
  n_cluster = 5L,
  knn = list()
) {
  checkmate::qassert(n_perm, "I1")
  checkmate::qassert(n_cluster, "I1")

  knn_params <- modifyList(
    params_knn_defaults(),
    modifyList(list(k = 15L, nn_max_iter = 15L), knn, keep.null = TRUE),
    keep.null = TRUE
  )

  c(
    list(
      n_perm = n_perm,
      n_cluster = n_cluster
    ),
    knn_params
  )
}

## hotspot ---------------------------------------------------------------------

#' Wrapper function for parameters for HotSpot
#'
#' @param model String. Model to use for modelling the GEX. One of
#' `c("danb", "bernoulli", "normal")`. Defaults to `"danb"`.
#' @param normalise Boolean. Shall the data be normalised. Defaults to `TRUE`.
#' @param knn List. Optional overrides for kNN parameters. See
#' [bixverse::params_knn_defaults()] for available parameters: `k`,
#' `knn_method`, `ann_dist`, `search_budget`, `n_trees`, `delta`,
#' `diversify_prob`, `ef_budget`, `m`, `ef_construction`, `ef_search`, `n_list`
#' and `n_probe`.
#'
#' @returns A list with the HotSpot parameters.
#'
#' @export
params_sc_hotspot <- function(
  model = c("danb", "normal", "bernoulli"),
  normalise = TRUE,
  knn = list(ann_dist = "cosine")
) {
  model <- match.arg(model)
  checkmate::qassert(normalise, "B1")

  knn_params <- modifyList(
    params_knn_defaults(),
    knn,
    keep.null = TRUE
  )

  list(
    knn_method = knn_params$knn_method,
    ann_dist = knn_params$ann_dist,
    k = knn_params$k,
    n_tree = knn_params$n_trees,
    search_budget = knn_params$search_budget,
    max_iter = knn_params$nn_max_iter,
    rho = knn_params$rho,
    delta = knn_params$delta,
    model = model,
    normalise = normalise
  )
}

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
#' `c("hnsw", "annoy", "nndescent", "ivf")`. Defaults to `"nndescent"`.
#' @param knn List. Optional overrides for kNN parameters. See
#' [bixverse::params_knn_defaults()] for available parameters: `k`,
#' `knn_method`, `ann_dist`, `search_budget`, `n_trees`, `delta`,
#' `diversify_prob`, `ef_budget`, `m`, `ef_construction`, `ef_search`, `n_list`
#' and `n_probe`. Note: `knn_method` cannot be `"exhaustive"` for MiloR as it
#' doesn't generate an index!
#'
#' @returns A list with the MiloR parameters.
#'
#' @export
params_sc_miloR <- function(
  prop = 0.2,
  k_refine = 20L,
  refinement_strategy = c("index", "approximate", "bruteforce"),
  index_type = c("nndescent", "ivf", "hnsw", "annoy"),
  knn = list()
) {
  refinement_strategy <- match.arg(refinement_strategy)
  index_type <- match.arg(index_type)
  checkmate::qassert(prop, "N1(0,1)")
  checkmate::qassert(k_refine, "I1")

  knn_params <- modifyList(
    params_knn_defaults(),
    knn,
    keep.null = TRUE
  )

  list(
    prop = prop,
    k_refine = k_refine,
    refinement_strategy = refinement_strategy,
    index_type = index_type,
    knn_method = knn_params$knn_method,
    ann_dist = knn_params$ann_dist,
    k = knn_params$k,
    n_trees = knn_params$n_trees,
    search_budget = knn_params$search_budget,
    nn_max_iter = knn_params$nn_max_iter,
    rho = knn_params$rho,
    delta = knn_params$delta
  )
}

## metacells -------------------------------------------------------------------

### meta cell (hdWGCNA) --------------------------------------------------------

#' Wrapper function for parameters for bootstrapped meta cell generation
#'
#' @description
#' This function generates parameters for the bootstrapped meta cell generation
#' based on hdWGCNA, see Morabito, et al., Cell Rep. Methods, 2023.
#'
#' @param max_shared Integer. Maximum number of allowed shared neighbours for
#' the meta cell to be considered. Defaults to `15L`.
#' @param target_no_metacells Integer. Target number of meta-cells to generate.
#' Defaults to `1000L`.
#' @param max_iter Integer. Maximum number of iterations for the algorithm.
#' Defaults to `5000L`.
#' @param knn List. Optional overrides for kNN parameters. See
#' [bixverse::params_knn_defaults()] for available parameters: `k`,
#' `knn_method`, `ann_dist`, `search_budget`, `n_trees`, `delta`,
#' `diversify_prob`, `ef_budget`, `m`, `ef_construction`, `ef_search`, `n_list`
#' and `n_probe`.
#'
#' @returns A list with the metacell parameters.
#'
#' @export
params_sc_bt_metacells <- function(
  max_shared = 15L,
  target_no_metacells = 1000L,
  max_iter = 5000L,
  knn = list()
) {
  checkmate::qassert(max_shared, "I1")
  checkmate::qassert(target_no_metacells, "I1")
  checkmate::qassert(max_iter, "I1")

  knn_params <- modifyList(
    params_knn_defaults(),
    modifyList(list(k = 25L, ann_dist = "cosine"), knn, keep.null = TRUE),
    keep.null = TRUE
  )

  c(
    list(
      max_shared = max_shared,
      target_no_metacells = target_no_metacells,
      max_iter = max_iter
    ),
    knn_params
  )
}

### sea cells ------------------------------------------------------------------

#' Wrapper function for the SEACells parameters
#'
#' @param n_sea_cells Integer. Number of SEA cells to detect.
#' @param max_fw_iters Integer. Maximum iterations for the Franke-Wolfe
#' algorithm. Defaults to `50L`.
#' @param convergence_epsilon Numeric. Convergence threshold. Algorithm stops
#' when RSS change < epsilon * RSS(0). Defaults to `1e-3`.
#' @param max_iter Integer. Maximum iterations to run SEACells for. Defaults to
#' `100L`.
#' @param min_iter Integer. Minimum iterations to run SEACells for. Defaults to
#' `10L`.
#' @param greedy_threshold Integer. Maximum number of cells before defaulting to
#' rapid random selection of archetypes. Defaults to `20000L`.
#' @param graph_building String. Graph building method. Defaults to `"union"`.
#' @param pruning Boolean. Shall tiny values be pruned during Franke-Wolfe
#' updates. This will reduce memory pressure and can be a good option on
#' large data sets. Defaults to `FALSE`.
#' @param pruning_threshold Float. If `pruning = TRUE` values below which
#' threshold shall be pruned.
#' @param n_landmarks Optional integer. If provided, it will use the Nystroem
#' extension during the archetype finding. Useful for larger data sets.
#' @param knn List. Optional overrides for kNN parameters. See
#' [bixverse::params_knn_defaults()] for available parameters: `k`,
#' `knn_method`, `ann_dist`, `search_budget`, `n_trees`, `delta`,
#' `diversify_prob`, `ef_budget`, `m`, `ef_construction`, `ef_search`, `n_list`
#' and `n_probe`.
#'
#' @returns A list with the SEACells parameters.
#'
#' @export
params_sc_seacells <- function(
  n_sea_cells,
  max_fw_iters = 50L,
  convergence_epsilon = 1e-3,
  max_iter = 100L,
  min_iter = 10L,
  greedy_threshold = 20000L,
  graph_building = "union",
  pruning = FALSE,
  pruning_threshold = 1e-7,
  n_landmarks = NULL,
  knn = list()
) {
  checkmate::qassert(n_sea_cells, "I1")
  checkmate::qassert(max_fw_iters, "I1")
  checkmate::qassert(convergence_epsilon, "N1")
  checkmate::qassert(max_iter, "I1")
  checkmate::qassert(min_iter, "I1")
  checkmate::qassert(greedy_threshold, "I1")
  checkmate::qassert(graph_building, "S1")
  checkmate::qassert(pruning, "B1")
  checkmate::qassert(pruning_threshold, "N1")
  checkmate::qassert(n_landmarks, c("0", "N1"))

  knn_params <- modifyList(
    params_knn_defaults(),
    modifyList(list(k = 25L), knn, keep.null = TRUE),
    keep.null = TRUE
  )

  c(
    list(
      n_sea_cells = n_sea_cells,
      max_fw_iters = max_fw_iters,
      convergence_epsilon = convergence_epsilon,
      max_iter = max_iter,
      min_iter = min_iter,
      greedy_threshold = greedy_threshold,
      graph_building = graph_building,
      pruning = pruning,
      pruning_threshold = pruning_threshold,
      n_landmarks = n_landmarks
    ),
    knn_params
  )
}

### supercell ------------------------------------------------------------------

#' Wrapper function for parameters for SuperCell generation
#'
#' @param walk_length Integer. Walk length for the Walktrap algorithm. Defaults
#' to `3L`.
#' @param graining_factor Numeric. Graining level of data (proportion of number
#' of single cells in the initial dataset to the number of metacells in the
#' final dataset). Defaults to `20.0`. (One meta cell per 20 cells.)
#' @param use_kernel Boolean. Shall a kernel function akin to MAGIC be applied
#' akin to the approach in SuperCell2, see Hérault, et al., bioRxiv, 2026 and
#' van Dijk, et al., Cell, 2018.
#' @param k_ith_neighbour Optional integer. The k-ith neighbour to use for
#' the kernel. Defaults to `k %/% 2`.
#' @param max_support Optional integer. Caps each cell's walk-probability vector
#' to its top entries by mass, bounding memory at ~`max_support * n_cells` on
#' large data. Makes the result an approximation. `NULL` (default) keeps the
#' walks exact.
#' @param knn List. Optional overrides for kNN parameters. See
#' [bixverse::params_knn_defaults()] for available parameters: `k`,
#' `knn_method`, `ann_dist`, `search_budget`, `n_trees`, `delta`,
#' `diversify_prob`, `ef_budget`, `m`, `ef_construction`, `ef_search`, `n_list`
#' and `n_probe`.
#'
#' @returns A list with the SuperCell parameters.
#'
#' @export
params_sc_supercell <- function(
  walk_length = 3L,
  graining_factor = 20.0,
  use_kernel = TRUE,
  k_ith = NULL,
  max_support = NULL,
  knn = list()
) {
  checkmate::qassert(walk_length, "I1")
  checkmate::qassert(graining_factor, "N1")
  checkmate::qassert(use_kernel, "B1")
  checkmate::qassert(k_ith, c("I1", "0"))

  knn_params <- modifyList(
    params_knn_defaults(),
    modifyList(list(k = 5L, ann_dist = "cosine"), knn, keep.null = TRUE),
    keep.null = TRUE
  )

  c(
    list(
      walk_length = walk_length,
      graining_factor = graining_factor,
      use_kernel = use_kernel,
      k_ith = k_ith,
      max_support = max_support
    ),
    knn_params
  )
}

## batch correction methods ----------------------------------------------------

### BBKNN ----------------------------------------------------------------------

#' Wrapper function for the BBKNN parameters
#'
#' @param neighbours_within_batch Integer. Number of neighbours to consider
#' per batch. Defaults to `3L`.
#' @param set_op_mix_ratio Numeric. Mixing ratio between union (1.0) and
#' intersection (0.0). Defaults to `1.0`.
#' @param local_connectivity Numeric. UMAP connectivity computation parameter,
#' how many nearest neighbours of each cell are assumed to be fully connected.
#' Defaults to `1.0`.
#' @param trim Optional integer. Trim the neighbours of each cell to these many
#' top connectivities. May help with population independence and improve the
#' tidiness of clustering. If `NULL`, it defaults to
#' `10 * neighbours_within_batch`.
#' @param knn List. Optional overrides for kNN parameters. See
#' [bixverse::params_knn_defaults()] for available parameters: `k`,
#' `knn_method`, `ann_dist`, `search_budget`, `n_trees`, `delta`,
#' `diversify_prob`, `ef_budget`, `m`, `ef_construction`, `ef_search`, `n_list`
#' and `n_probe`.
#'
#' @returns A list with the BBKNN parameters.
#'
#' @export
params_sc_bbknn <- function(
  neighbours_within_batch = 3L,
  set_op_mix_ratio = 1.0,
  local_connectivity = 1.0,
  trim = NULL,
  knn = list()
) {
  checkmate::qassert(neighbours_within_batch, "I1")
  checkmate::qassert(set_op_mix_ratio, "N1[0,1]")
  checkmate::qassert(local_connectivity, "N1")
  checkmate::qassert(trim, c("0", "I1"))

  knn_params <- modifyList(
    params_knn_defaults(),
    modifyList(list(k = neighbours_within_batch * 2L), knn, keep.null = TRUE),
    keep.null = TRUE
  )

  c(
    list(
      neighbours_within_batch = neighbours_within_batch,
      set_op_mix_ratio = set_op_mix_ratio,
      local_connectivity = local_connectivity,
      trim = trim
    ),
    knn_params
  )
}

### fastMNN --------------------------------------------------------------------

#' Wrapper function for the fastMNN parameters
#'
#' @param ndist Numeric. Number of median distances for the tricube kernel
#' bandwidth. Defaults to `3.0`.
#' @param cos_norm Logical. Apply cosine normalisation before computing
#' distances. Defaults to `TRUE`.
#' @param no_pcs Integer. Number of PCs to use for MNN calculations.
#' Defaults to `30L`.
#' @param sparse_svd Boolean. Shall the sparse SVD be used.
#' @param knn List. Optional overrides for kNN parameters. See
#' [bixverse::params_knn_defaults()] for available parameters: `k`,
#' `knn_method`, `ann_dist`, `search_budget`, `n_trees`, `delta`,
#' `diversify_prob`, `ef_budget`, `m`, `ef_construction`, `ef_search`, `n_list`
#' and `n_probe`.
#' @param pca Named list. Parameters to feed through to the optional
#' recalculation of the PCA, see [params_sc_pca()].
#'
#' @returns A list with the fastMNN parameters.
#'
#' @export
params_sc_fastmnn <- function(
  ndist = 3.0,
  cos_norm = TRUE,
  no_pcs = 30L,
  sparse_svd = TRUE,
  knn = list(k = 20L),
  pca = params_sc_pca()
) {
  checkmate::qassert(ndist, "N1(0,)")
  checkmate::qassert(cos_norm, "B1")
  checkmate::qassert(no_pcs, "I1")
  checkmate::qassert(random_svd, "B1")
  checkmate::qassert(sparse_svd, "B1")

  knn_params <- modifyList(
    params_knn_defaults(),
    modifyList(list(k = 20L, ann_dist = "cosine"), knn, keep.null = TRUE),
    keep.null = TRUE
  )

  c(
    list(
      ndist = ndist,
      cos_norm = cos_norm,
      no_pcs = no_pcs,
      sparse_svd = sparse_svd
    ),
    knn_params,
    pca
  )
}

### harmony --------------------------------------------------------------------

#' Default parameters for Harmony batch correction
#'
#' @param k Optional integer. Number of clusters for k-means clustering. If
#' not provided, it will be automatically determined as
#' `min(round(N / 30), 100)`.
#' @param sigma Numeric vector. Per-cluster diversity weights. Either a single
#' value (broadcast to all clusters) or a vector of length k.
#' @param theta Numeric vector. Per-variable diversity penalties. Either a
#' single value (broadcast to all variables) or a vector of length equal to the
#' number of batch variables.
#' @param lambda Numeric vector. Ridge regression penalty for the linear model.
#' Typically a single value that is broadcast to all design matrix columns.
#' @param block_size Numeric. Fraction of cells to update per block during
#' optimisation (0.0-1.0). Lower values reduce memory usage but increase
#' computation time.
#' @param max_iter_kmeans Integer. Maximum number of k-means iterations per
#' Harmony round.
#' @param max_iter_harmony Integer. Maximum number of Harmony outer iterations.
#' @param epsilon_kmeans Numeric. Convergence threshold for k-means clustering.
#' Stops when the relative change in cluster assignments falls below this value.
#' @param epsilon_harmony Numeric. Convergence threshold for Harmony. Stops when
#' the relative change in the objective function falls below this value.
#' @param window_size Integer. Number of previous iterations to consider when
#' checking convergence.
#' @param kmeans List. Optional overrides for the k-means clustering algorithm
#' Possible parameters are `"k_means_iter"`, `"k_means_init"`, `"gemm"` and
#' `"hamerly"`, see [params_kmeans_defaults()].
#'
#' @return A list with the parameters.
#'
#' @export
params_sc_harmony <- function(
  k = NULL,
  sigma = 0.1,
  theta = 2.0,
  lambda = 1.0,
  block_size = 0.2,
  max_iter_kmeans = 20L,
  max_iter_harmony = 10L,
  epsilon_kmeans = 1e-5,
  epsilon_harmony = 1e-4,
  window_size = 2L,
  kmeans = list()
) {
  # checks
  checkmate::qassert(k, c("I1[1,)", "0"))
  checkmate::qassert(sigma, "N+[0,)")
  checkmate::qassert(theta, "N+[0,)")
  checkmate::qassert(lambda, "N+[0,)")
  checkmate::qassert(block_size, "N1(0,1]")
  checkmate::qassert(max_iter_kmeans, "I1[1,)")
  checkmate::qassert(max_iter_harmony, "I1[1,)")
  checkmate::qassert(epsilon_kmeans, "N1(0,)")
  checkmate::qassert(epsilon_harmony, "N1(0,)")
  checkmate::qassert(window_size, "I1[1,)")

  kmeans_params <- modifyList(
    params_kmeans_defaults(),
    kmeans,
    keep.null = TRUE
  )

  c(
    list(
      k = k,
      sigma = sigma,
      theta = theta,
      lambda = lambda,
      block_size = block_size,
      max_iter_kmeans = max_iter_kmeans,
      max_iter_harmony = max_iter_harmony,
      epsilon_kmeans = epsilon_kmeans,
      epsilon_harmony = epsilon_harmony,
      window_size = window_size
    ),
    kmeans_params
  )
}

### harmony (version 2) --------------------------------------------------------

#' Default parameters for Harmony v2 batch correction
#'
#' @param k Optional integer. Number of clusters for k-means clustering. If
#' not provided, it will be automatically determined as
#' `min(round(N / 30), 100)`.
#' @param sigma Numeric vector. Per-cluster diversity weights. Either a single
#' value (broadcast to all clusters) or a vector of length k.
#' @param theta Numeric vector. Per-variable diversity penalties. Either a
#' single value (broadcast to all variables) or a vector of length equal to the
#' number of batch variables.
#' @param lambda Numeric vector. Ridge regression penalty for the linear model.
#' Typically a single value that is broadcast to all design matrix columns.
#' Ignored when `use_dynamic_lambda = TRUE`.
#' @param block_size Numeric. Fraction of cells to update per block during
#' optimisation (0.0-1.0). Lower values reduce memory usage but increase
#' computation time.
#' @param max_iter_kmeans Integer. Maximum number of k-means iterations per
#' Harmony round.
#' @param max_iter_harmony Integer. Maximum number of Harmony outer iterations.
#' @param epsilon_kmeans Numeric. Convergence threshold for k-means clustering.
#' Stops when the relative change in cluster assignments falls below this value.
#' @param epsilon_harmony Numeric. Convergence threshold for Harmony. Stops when
#' the relative change in the objective function falls below this value.
#' @param window_size Integer. Number of previous iterations to consider when
#' checking convergence.
#' @param alpha Numeric. Scaling factor for dynamic lambda estimation. Must be
#' in (0, 1). Only relevant when `use_dynamic_lambda = TRUE`.
#' @param tau Numeric. Scaling factor for theta based on batch size. A value of
#' 0 disables batch-size scaling of theta.
#' @param batch_proportion_cutoff Numeric. Cutoff for pruning batches with small
#' proportions during ridge regression.
#' @param use_dynamic_lambda Boolean. If `TRUE`, lambda is estimated dynamically
#' per cluster instead of using the fixed `lambda` value.
#' @param kmeans List. Optional overrides for the k-means clustering algorithm
#' Possible parameters are `"k_means_iter"`, `"k_means_init"`, `"gemm"` and
#' `"hamerly"`, see [params_kmeans_defaults()].
#'
#' @return A list with the parameters.
#'
#' @export
params_sc_harmony_v2 <- function(
  k = NULL,
  sigma = 0.1,
  theta = 2.0,
  lambda = 1.0,
  block_size = 0.2,
  max_iter_kmeans = 4L,
  max_iter_harmony = 10L,
  epsilon_kmeans = 1e-3,
  epsilon_harmony = 1e-2,
  window_size = 3L,
  alpha = 0.2,
  tau = 0.0,
  batch_proportion_cutoff = 1e-5,
  use_dynamic_lambda = FALSE,
  kmeans = list()
) {
  # checks
  checkmate::qassert(k, c("I1[1,)", "0"))
  checkmate::qassert(sigma, "N+[0,)")
  checkmate::qassert(theta, "N+[0,)")
  checkmate::qassert(lambda, "N+[0,)")
  checkmate::qassert(block_size, "N1(0,1]")
  checkmate::qassert(max_iter_kmeans, "I1[1,)")
  checkmate::qassert(max_iter_harmony, "I1[1,)")
  checkmate::qassert(epsilon_kmeans, "N1(0,)")
  checkmate::qassert(epsilon_harmony, "N1(0,)")
  checkmate::qassert(window_size, "I1[1,)")
  checkmate::qassert(alpha, "N1(0,1)")
  checkmate::qassert(tau, "N1[0,)")
  checkmate::qassert(batch_proportion_cutoff, "N1(0,)")
  checkmate::qassert(use_dynamic_lambda, "B1")

  kmeans_params <- modifyList(
    params_kmeans_defaults(),
    kmeans,
    keep.null = TRUE
  )

  c(
    list(
      k = k,
      sigma = sigma,
      theta = theta,
      lambda = lambda,
      block_size = block_size,
      max_iter_kmeans = max_iter_kmeans,
      max_iter_harmony = max_iter_harmony,
      epsilon_kmeans = epsilon_kmeans,
      epsilon_harmony = epsilon_harmony,
      window_size = window_size,
      alpha = alpha,
      tau = tau,
      batch_proportion_cutoff = batch_proportion_cutoff,
      use_dynamic_lambda = use_dynamic_lambda
    ),
    kmeans_params
  )
}

## scenic ----------------------------------------------------------------------

### regression learner params --------------------------------------------------

#' Default parameters for the SCENIC RandomForest regression learner
#'
#' @return A list with the following parameters:
#' \itemize{
#'  \item n_trees - Integer. Number of trees to build. Defaults to `250L`.
#'  \item min_samples_leaf - Integer. Minimum number of samples required at a
#'  leaf node. Defaults to `50L`.
#'  \item n_features_split - Integer. Number of features considered per split.
#'  `0L` resolves to `sqrt(n_features)` at runtime. Defaults to `0L`.
#'  \item subsample_rate - Numeric. Fraction of samples to draw per tree.
#'  Defaults to `0.632`.
#'  \item bootstrap - Logical. Whether to sample with replacement. Defaults to
#'  `FALSE`.
#'  \item max_depth - Integer. Maximum depth of each tree. Defaults to `8L`.
#'  \item subsample_frac - Optional numeric. Fraction of cells to subsample per
#'  tree. If set, overrides `subsample_rate`. Defaults to `NULL`.
#' }
#'
#' @export
params_scenic_random_forest_defaults <- function() {
  list(
    n_trees = 250L,
    min_samples_leaf = 50L,
    n_features_split = 0L,
    subsample_rate = 0.632,
    bootstrap = FALSE,
    max_depth = 8L,
    subsample_frac = NULL
  )
}

#' Default parameters for the SCENIC ExtraTrees regression learner
#'
#' @return A list with the following parameters:
#' \itemize{
#'  \item n_trees - Integer. Number of trees to build. Defaults to `500L`.
#'  \item min_samples_leaf - Integer. Minimum number of samples required at a
#'  leaf node. Defaults to `50L`.
#'  \item n_features_split - Integer. Number of features considered per split.
#'  `0L` resolves to `sqrt(n_features)` at runtime. Defaults to `0L`.
#'  \item n_thresholds - Integer. Number of random thresholds to evaluate per
#'  feature per node. Defaults to `1L`.
#'  \item max_depth - Integer. Maximum depth of each tree. Defaults to `8L`.
#'  \item subsample_frac - Optional numeric. Fraction of cells to subsample per
#'  tree. Defaults to `NULL`.
#' }
#'
#' @export
params_scenic_extra_trees_defaults <- function() {
  list(
    n_trees = 500L,
    min_samples_leaf = 50L,
    n_features_split = 0L,
    n_thresholds = 1L,
    max_depth = 8L,
    subsample_frac = NULL
  )
}

#' Default parameters for the SCENIC GradientBoosting (GRNBoost2) regression
#' learner
#'
#' @return A list with the following parameters:
#' \itemize{
#'  \item n_trees_max - Integer. Maximum number of boosting rounds. Early
#'  stopping usually triggers well before this limit. Defaults to `1000L`.
#'  \item learning_rate - Numeric. Shrinkage applied to each tree's
#'  predictions. Defaults to `0.01`.
#'  \item max_depth - Integer. Maximum depth of each tree. Shallow trees
#'  (3-5) work best for GBM. Defaults to `3L`.
#'  \item min_samples_leaf - Integer. Minimum number of training samples
#'  required at a leaf node. Defaults to `50L`.
#'  \item early_stop_window - Integer. Number of recent OOB improvements to
#'  average for the early stopping criterion. Stops when the rolling average
#'  drops to zero or below. Defaults to `25L`.
#'  \item subsample_rate - Numeric. Fraction of samples used for training
#'  each tree. The complement forms the OOB set. Defaults to `0.9`.
#'  \item n_features_split - Integer. Number of features to evaluate per
#'  split. `0L` means all features (recommended with histogram subtraction).
#'  Defaults to `0L`.
#' }
#'
#' @export
params_scenic_gradient_boosting_defaults <- function() {
  list(
    n_trees_max = 1000L,
    learning_rate = 0.01,
    max_depth = 3L,
    min_samples_leaf = 50L,
    early_stop_window = 25L,
    subsample_rate = 0.9,
    n_features_split = 0L
  )
}

### main params ----------------------------------------------------------------

#' Constructor for SCENIC parameters
#'
#' @param min_counts Integer. Minimum total counts a gene needs to be included
#' in the analysis. Defaults to `50L`.
#' @param min_cells Numeric. Minimum proportion of cells (between 0 and 1) that
#' must express a gene for it to be considered. Defaults to `0.03`.
#' @param learner_type Character. Regression learner to use. One of
#' `"randomforest"`, `"extratrees"`, or `"grnboost2"`. Defaults to
#' `"randomforest"`.
#' @param gene_batch_strategy Character. Strategy for grouping target genes into
#' batches. One of `"random"` or `"correlated"`. Only used for `"randomforest"`
#' and `"extratrees"` learners; ignored for `"grnboost2"`. Defaults to
#' `"correlated"`.
#' @param gene_batch_size Optional integer. Number of genes per batch. If `NULL`
#' (default), the batch size is determined automatically. Ignored for
#' `"grnboost2"`.
#' @param n_pcs Integer. Number of PCs to use for the correlated gene batch
#' strategy. Defaults to `50L`.
#' @param n_subsample Integer. Cell subsampling threshold for the correlated
#' gene batch strategy. If the number of cells meets or exceeds this value,
#' `n_subsample` cells are randomly selected prior to running randomised SVD.
#' Defaults to `100000L`.
#' @param learner_params List. Optional overrides for the regression learner
#' parameters. For `"randomforest"`, see
#' [bixverse::params_scenic_random_forest_defaults()]. For `"extratrees"`, see
#' [bixverse::params_scenic_extra_trees_defaults()]. For `"grnboost2"`, see
#' [bixverse::params_scenic_gradient_boosting_defaults()].
#'
#' @returns A named flat list with all SCENIC parameters.
#'
#' @export
params_scenic <- function(
  min_counts = 50L,
  min_cells = 0.03,
  learner_type = "randomforest",
  gene_batch_strategy = "correlated",
  gene_batch_size = NULL,
  n_pcs = 50L,
  n_subsample = 100000L,
  learner_params = list()
) {
  checkmate::qassert(min_counts, "I1[1,)")
  checkmate::qassert(min_cells, "N1(0,1]")
  checkmate::assert_choice(
    learner_type,
    c("randomforest", "extratrees", "grnboost2")
  )
  checkmate::assert_choice(gene_batch_strategy, c("random", "correlated"))
  if (!is.null(gene_batch_size)) {
    checkmate::qassert(gene_batch_size, "I1[1,)")
  }
  checkmate::qassert(n_pcs, "I1[1,)")
  checkmate::qassert(n_subsample, "I1[1,)")

  learner_defaults <- switch(
    learner_type,
    extratrees = params_scenic_extra_trees_defaults(),
    grnboost2 = params_scenic_gradient_boosting_defaults(),
    params_scenic_random_forest_defaults()
  )

  params <- c(
    list(
      min_counts = min_counts,
      min_cells = min_cells,
      learner_type = learner_type,
      gene_batch_strategy = gene_batch_strategy,
      gene_batch_size = gene_batch_size,
      n_pcs = n_pcs,
      n_subsample = n_subsample
    ),
    modifyList(learner_defaults, learner_params, keep.null = TRUE)
  )

  params
}

## meld ------------------------------------------------------------------------

#' Constructor for MELD parameters
#'
#' @param beta Numeric. Smoothing strength; larger values produce smoother
#' densities. Must be strictly positive. Defaults to `60.0`.
#' @param offset Numeric. Shift of the filter centre in the rescaled spectrum.
#' Must be in `[0, 1]`. Defaults to `0.0`.
#' @param order Numeric. Filter falloff sharpness; larger values approach a
#' square low-pass. Must be strictly positive. Defaults to `1.0`.
#' @param filter Character. Filter family to use. One of `"heat"` or
#' `"laplacian"`. Defaults to `"heat"`.
#' @param chebyshev_order Integer. Number of Chebyshev coefficients (polynomial
#' terms). Must be >= 2. Defaults to `50L`.
#' @param lap_type Character. Type of Laplacian to use for spectral filtering.
#' One of `"combinatorial"` or `"normalised"`. Defaults to `"combinatorial"`.
#' @param normalise_indicators Logical. If `TRUE`, each column of the indicator
#' matrix is divided by its column sum before filtering, making
#' cross-condition densities comparable regardless of cells-per-condition.
#' Defaults to `TRUE`.
#' @param knn List. Optional overrides for kNN parameters. See
#' [bixverse::params_knn_defaults()] for available parameters: `k`,
#' `knn_method`, `ann_dist`, `search_budget`, `n_trees`, `delta`,
#' `diversify_prob`, `ef_budget`, `m`, `ef_construction`, `ef_search`,
#' `n_list` and `n_probe`.
#'
#' @returns A named flat list with all MELD parameters.
#'
#' @export
params_meld <- function(
  beta = 60.0,
  offset = 0.0,
  order = 1.0,
  filter = "heat",
  chebyshev_order = 50L,
  lap_type = "combinatorial",
  normalise_indicators = TRUE,
  knn = list()
) {
  checkmate::qassert(beta, "N1(0,)")
  checkmate::qassert(offset, "N1[0,1]")
  checkmate::qassert(order, "N1(0,)")
  checkmate::assert_choice(filter, c("heat", "laplacian"))
  checkmate::qassert(chebyshev_order, "I1[2,)")
  checkmate::assert_choice(lap_type, c("combinatorial", "normalised"))
  checkmate::qassert(normalise_indicators, "B1")

  params <- list(
    knn = modifyList(params_knn_defaults(), knn, keep.null = TRUE),
    beta = beta,
    offset = offset,
    order = order,
    filter = filter,
    chebyshev_order = chebyshev_order,
    lap_type = lap_type,
    normalise_indicators = normalise_indicators
  )

  params <- purrr::list_flatten(params, name_spec = "{inner}")

  params
}

## multi-modal -----------------------------------------------------------------

### wnn ------------------------------------------------------------------------

#' Wrapper function for WNN parameters
#'
#' @param k_nn Integer. Final number of multimodal neighbours per cell. Defaults
#' to `20L`.
#' @param knn_range Integer. Candidate pool size per modality. Each cell's kNN
#' input must contain at least this many neighbours. Defaults to `100L`.
#' @param sigma_method String. Bandwidth method. One of
#' `c("snn_farthest", "sigma_idx")`. Defaults to `"snn_farthest"`.
#' @param sigma_idx Integer. `"sigma_idx"` only: 0-based kNN index for
#' bandwidth. Defaults to `19L` (i.e. `k_nn - 1`).
#' @param snn_type String. sNN type. One of `c("full_connection", "limited")`.
#' The limited version only considers edges that exist in the kNN. Defaults to
#' `"full_connection"`.
#' @param s_nn Integer. `"snn_farthest"` only: kNN size used to build the SNN
#' graph. Defaults to `20L`.
#' @param sd_scale Numeric. Multiplier on sigma. Defaults to `1.0`.
#' @param kernel_power Numeric. Kernel exponent power. Defaults to `1.0`.
#' @param cross_const Numeric. Cross-modality kernel stabiliser. Defaults to
#' `1e-4`.
#' @param sigma_floor Numeric. Minimum sigma value (avoids division by zero).
#' Defaults to `1e-8`.
#' @param knn List. Optional overrides for kNN parameters. See
#' [bixverse::params_knn_defaults()] for available parameters: `k`,
#' `knn_method`, `ann_dist`, `search_budget`, `n_trees`, `delta`,
#' `diversify_prob`, `ef_budget`, `m`, `ef_construction`, `ef_search`, `n_list`
#' and `n_probe`.
#'
#' @returns A list with the WNN parameters.
#'
#' @export
params_sc_wnn <- function(
  k_nn = 20L,
  knn_range = 200L,
  sigma_method = c("snn_farthest", "sigma_idx"),
  sigma_idx = 19L,
  snn_type = c("full_connection", "limited"),
  s_nn = 20L,
  sd_scale = 1.0,
  kernel_power = 1.0,
  cross_const = 1e-4,
  sigma_floor = 1e-8,
  knn = list()
) {
  sigma_method <- match.arg(sigma_method)
  snn_type <- match.arg(snn_type)

  checkmate::qassert(k_nn, "I1[1,)")
  checkmate::qassert(knn_range, "I1[1,)")
  checkmate::assertChoice(sigma_method, c("snn_farthest", "sigma_idx"))
  checkmate::qassert(sigma_idx, "I1[0,)")
  checkmate::assertChoice(snn_type, c("full_connection", "limited"))
  checkmate::qassert(s_nn, "I1[1,)")
  checkmate::qassert(sd_scale, "N1(0,)")
  checkmate::qassert(kernel_power, "N1(0,)")
  checkmate::qassert(cross_const, "N1[0,)")
  checkmate::qassert(sigma_floor, "N1(0,)")

  knn_params <- modifyList(
    params_knn_defaults(),
    knn,
    keep.null = TRUE
  )

  c(
    list(
      k_nn = k_nn,
      knn_range = knn_range,
      sigma_method = sigma_method,
      sigma_idx = sigma_idx,
      snn_type = snn_type,
      s_nn = s_nn,
      sd_scale = sd_scale,
      kernel_power = kernel_power,
      cross_const = cross_const,
      sigma_floor = sigma_floor
    ),
    knn_params
  )
}
