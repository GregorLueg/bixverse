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
#' @param gbm_n_trees Integer. Maximum number of boosting rounds for the GBM
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
#' @param cxds_genes Optional integer. Number of CXDS genes to consider. If not
#' provided, defaults to `500L`.
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
  gbm_n_trees = 200L,
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
  checkmate::qassert(gbm_n_trees, "I1[1,)")
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
    gbm_n_trees = gbm_n_trees,
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
#' Track this against `k` rather than leaving it: the threshold is a share of
#' the neighbourhood, so the same value prunes far harder at a larger `k`.
#' Over-pruning fails quietly, in that you still get a clustering, but cells
#' left with too few shared neighbours drop out as singleton communities, which
#' then show up downstream as one-cell clusters with inflated
#' [bixverse::run_paga_sc()] connectivities.
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
  knn = list()
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

## aucell ----------------------------------------------------------------------

#' Wrapper function for parameters for AUCell
#'
#' @description
#' The three statistics consume the same within-cell ranking, but weight it
#' very differently. `"wilcox"` is a pure function of the gene set's rank sum,
#' so a gene at rank 2 and one at rank 200 count for almost the same thing.
#' `"recovery"` and `"ap"` are top-heavy.
#'
#' For SCENIC use `"recovery"`, the default. `"wilcox"` is a bixverse addition
#' and its flatter score does not separate into on/off populations, so it
#' binarises badly.
#'
#' Ranking ties are averaged (midranks) rather than broken at random the way
#' AUCell does it, which is why there is no need for SCENIC's trick of setting
#' the cutoff to the 1st percentile of genes detected per cell. Undetected
#' genes all collapse onto one rank well outside any sensible `max_rank`.
#'
#' Note the recovery AUC is normalised by `max_rank * length(gene_set)`,
#' matching pySCENIC and AUCell's `normAUC = FALSE`. Modern AUCell divides by
#' the attainable maximum instead, so absolute values differ by a per-gene-set
#' constant. Cell ordering within a regulon is unaffected.
#'
#' @param auc_type String. Which statistic to calculate. One of
#' `c("wilcox", "recovery", "ap")`. `"wilcox"` is the AUC derived from the
#' Mann-Whitney U statistic over the full ranking, with the null at 0.5 for any
#' gene set size. `"recovery"` is the recovery-curve AUC under a rank cutoff,
#' i.e. the actual AUCell statistic of Aibar, et al. `"ap"` is average
#' precision, the most top-heavy of the three, but its null tracks the gene set
#' prevalence so raw values are not comparable across gene sets of different
#' size unless `standardise` is on. Defaults to `"recovery"`.
#' @param max_rank Optional numeric. Rank cutoff for `"recovery"`, counted from
#' the top of the within-cell ranking. If `NULL`, resolves to the top 5% of the
#' gene universe, following Aibar, et al. Ignored by the other two statistics.
#' @param standardise Boolean. Shall each gene set's scores be z-scored across
#' the cells. This is what makes `"ap"` comparable across gene sets of
#' different size. Defaults to `FALSE`.
#'
#' @returns A list with the AUCell parameters.
#'
#' @references Aibar, et al., Nat Methods, 2017
#'
#' @export
params_sc_aucell <- function(
  auc_type = c("recovery", "wilcox", "ap"),
  max_rank = NULL,
  standardise = FALSE
) {
  auc_type <- match.arg(auc_type)
  checkmate::assertChoice(auc_type, c("recovery", "wilcox", "ap"))
  checkmate::qassert(max_rank, c("N1[1,)", "0"))
  checkmate::qassert(standardise, "B1")

  # Rust parses this one with as_real(), so an integer would silently be
  # dropped and fall back to the automatic cutoff
  list(
    auc_type = auc_type,
    max_rank = if (is.null(max_rank)) NULL else as.double(max_rank),
    standardise = standardise
  )
}

## scenic binarisation ---------------------------------------------------------

#' Wrapper function for parameters for the SCENIC binarisation
#'
#' @description
#' Each regulon gets its own threshold. A two-component Gaussian mixture is
#' fitted and compared against a single Gaussian by BIC; if the mixture wins,
#' the threshold is the kernel density minimum between the two component means,
#' otherwise it falls back to `mean + 2 * sd`. This follows pySCENIC. AUCell
#' fits six candidates and then lets the density trough override all of them
#' whenever one exists, so the two land in much the same place.
#'
#' Turn `bw_adjust` up if shallow wobbles in the density are being picked up as
#' troughs. AUCell effectively runs at `2`.
#'
#' @param bw_adjust Float. Multiplier on the Silverman bandwidth of the kernel
#' density estimate. Higher values smooth more. Defaults to `1`.
#' @param n_grid Integer. Number of points at which the density is evaluated
#' between the two component means. Defaults to `512L`.
#' @param n_bins Integer. Number of histogram bins used to approximate the
#' density. Defaults to `512L`.
#'
#' @returns A list with the binarisation parameters.
#'
#' @references Aibar, et al., Nat Methods, 2017
#'
#' @export
params_scenic_binarise <- function(
  bw_adjust = 1,
  n_grid = 512L,
  n_bins = 512L
) {
  checkmate::qassert(bw_adjust, "N1(0,)")
  checkmate::qassert(n_grid, "I1[3,)")
  checkmate::qassert(n_bins, "I1[2,)")

  # Rust parses these with as_real(), so integers need to go over as doubles
  list(
    bw_adjust = as.double(bw_adjust),
    n_grid = as.double(n_grid),
    n_bins = as.double(n_bins)
  )
}

## hotspot ---------------------------------------------------------------------

#' Wrapper function for parameters for HotSpot
#'
#' @description
#' `weighted_graph` controls how the kNN distances become edge weights. The
#' default of `FALSE` follows the reference implementation: the distances only
#' decide who is a neighbour and every retained edge weighs one. Set it to
#' `TRUE` for the Gaussian kernel, whose width is the
#' `ceil(k / neighborhood_factor)`-th neighbour distance.
#'
#' Whether the distances need squaring is derived from the metric, so it is not
#' a parameter here. When a pre-computed kNN graph is handed to the method, the
#' metric stored on that graph wins over `ann_dist`.
#'
#' @param model String. Model to use for modelling the GEX. One of
#' `c("danb", "bernoulli", "normal")`. Defaults to `"danb"`.
#' @param normalise Boolean. Shall the data be normalised. Defaults to `TRUE`.
#' @param weighted_graph Boolean. Shall the Gaussian kernel be applied to the
#' neighbour distances. Defaults to `FALSE`.
#' @param neighborhood_factor Float. Kernel width is the
#' `ceil(k / neighborhood_factor)`-th neighbour distance. Only read when
#' `weighted_graph = TRUE`. Defaults to `3`.
#' @param knn List. Optional overrides for kNN parameters. See
#' [bixverse::params_knn_defaults()] for available parameters: `k`,
#' `knn_method`, `ann_dist`, `search_budget`, `n_trees`, `delta`,
#' `diversify_prob`, `ef_budget`, `m`, `ef_construction`, `ef_search`, `n_list`
#' and `n_probe`.
#'
#' @returns A list with the HotSpot parameters.
#'
#' @references DeTomaso and Yosef, Cell Systems, 2021
#'
#' @export
params_sc_hotspot <- function(
  model = c("danb", "normal", "bernoulli"),
  normalise = TRUE,
  weighted_graph = FALSE,
  neighborhood_factor = 3.0,
  knn = list()
) {
  model <- match.arg(model)
  checkmate::assertChoice(model, c("danb", "normal", "bernoulli"))
  checkmate::qassert(normalise, "B1")
  checkmate::qassert(weighted_graph, "B1")
  checkmate::qassert(neighborhood_factor, "R1(0,)")

  knn_params <- modifyList(
    params_knn_defaults(),
    knn,
    keep.null = TRUE
  )

  c(
    list(
      model = model,
      normalise = normalise,
      weighted_graph = weighted_graph,
      neighborhood_factor = neighborhood_factor
    ),
    knn_params
  )
}

## dialogue --------------------------------------------------------------------

#' Wrapper function for the DIALOGUE decomposition parameters
#'
#' @description
#' Stage one of DIALOGUE: the penalised matrix decomposition that turns the
#' per-cell-type features into multicellular programmes, and the provisional
#' gene signatures that come off it.
#'
#' @details
#' The defaults follow upstream's `DLG.get.param`. Two knobs are worth thinking
#' about before anything else. `k` is how many programmes you are asking for,
#' and there is no sweep to help you pick it. `n_permutations` sets the
#' resolution of the empirical p-value: with the default of `100` the smallest
#' p you can observe is `0.01`, so lower it for a quick look and leave it alone
#' for anything you intend to believe.
#'
#' `averaging` is exposed and honoured here. Upstream takes the same argument
#' and then ignores it, hard-coding column medians, so `"median"` is what every
#' published DIALOGUE run actually used.
#'
#' @param k Integer. Number of multicellular programmes to extract. Must be at
#' least 1.
#' @param n_permutations Integer. Permutations backing the empirical p-value
#' per programme. Must be at least 2.
#' @param extra_sparse Boolean. Tune the L1 bound by permutation instead of
#' fixing it at `sqrt(p_1) / 2`. Costs ten more fits per permutation.
#' @param abn_c Integer. Minimum cells a sample must contribute, within a cell
#' type, before it counts towards the feature-level ANOVA.
#' @param p_anova Numeric. BH-adjusted ANOVA cutoff for keeping a feature. Must
#' be in `(0, 1]`.
#' @param centre Boolean. Centre and scale the sample-level feature matrix,
#' then winsorise it.
#' @param cap Numeric. Winsorising tail fraction applied to each column. Must
#' be in `[0, 0.5)`.
#' @param spatial Boolean. Spatial data: skip the ANOVA feature filter
#' entirely. Niches are small, so a feature need not vary across them to be
#' real.
#' @param n_genes Integer. Genes taken per programme per direction when
#' building a signature.
#' @param min_ci Numeric. Minimum absolute correlation for a gene to enter a
#' signature. Must be in `[0, 1]`.
#' @param averaging String. One of `c("median", "mean")`. How cell-level
#' features are collapsed per sample.
#' @param mcp_assignment_p Numeric. Empirical p below which a cell type pair
#' counts as connected when deciding which cell types a programme spans. Must
#' be in `(0, 1]`.
#' @param seed Integer. Seed for the permutation null.
#'
#' @returns A list with the stage one DIALOGUE parameters.
#'
#' @references Jerby-Arnon & Regev, Nature Biotechnology, 2022
#'
#' @export
params_dialogue_pmd <- function(
  k = 2L,
  n_permutations = 100L,
  extra_sparse = FALSE,
  abn_c = 15L,
  p_anova = 0.05,
  centre = TRUE,
  cap = 0.01,
  spatial = FALSE,
  n_genes = 200L,
  min_ci = 0.05,
  averaging = c("median", "mean"),
  mcp_assignment_p = 0.1,
  seed = 1234L
) {
  averaging <- match.arg(averaging)

  # checks
  checkmate::qassert(k, "I1[1,)")
  checkmate::qassert(n_permutations, "I1[2,)")
  checkmate::qassert(extra_sparse, "B1")
  checkmate::qassert(abn_c, "I1[0,)")
  checkmate::qassert(p_anova, "N1(0,1]")
  checkmate::qassert(centre, "B1")
  checkmate::qassert(cap, "N1[0,0.5)")
  checkmate::qassert(spatial, "B1")
  checkmate::qassert(n_genes, "I1[1,)")
  checkmate::qassert(min_ci, "N1[0,1]")
  checkmate::assertChoice(averaging, c("median", "mean"))
  checkmate::qassert(mcp_assignment_p, "N1(0,1]")
  checkmate::qassert(seed, "I1")

  list(
    k = k,
    n_permutations = n_permutations,
    extra_sparse = extra_sparse,
    abn_c = abn_c,
    p_anova = p_anova,
    centre = centre,
    cap = cap,
    spatial = spatial,
    n_genes = n_genes,
    min_ci = min_ci,
    averaging = averaging,
    mcp_assignment_p = mcp_assignment_p,
    seed = seed
  )
}

#' Wrapper function for the DIALOGUE mixed model parameters
#'
#' @description
#' Stage two of DIALOGUE: for every ordered pair of cell types and every
#' candidate gene, a random-intercept mixed model over samples asking whether a
#' cell's own programme score tracks the partner cell type's expression of that
#' gene in the same sample.
#'
#' @details
#' This stage dominates the runtime, and `satterthwaite` is the knob that
#' decides how badly. Turning it off falls back to the residual count for the
#' denominator degrees of freedom, which is far cheaper and barely differs once
#' a cell type has thousands of cells.
#'
#' `use_cell_quality` conditions on the cell's own quality covariate, which
#' stage one has already regressed out of the scores by ordinary least squares.
#' The default conditions on it twice, because upstream does.
#'
#' @param min_cells_per_sample Integer. Minimum cells a sample must contribute,
#' in *both* cell types of a pair, before it takes part in that pair's models.
#' @param use_tme_qc Boolean. Include the partner cell type's mean quality in
#' that sample as a fixed effect. Upstream's `tme.qc`.
#' @param use_cell_quality Boolean. Include the responding cell's own quality
#' as a fixed effect. Upstream's `cellQ`.
#' @param satterthwaite Boolean. Compute Satterthwaite denominator degrees of
#' freedom, as `lmerTest` does.
#'
#' @returns A list with the stage two DIALOGUE parameters.
#'
#' @references Jerby-Arnon & Regev, Nature Biotechnology, 2022
#'
#' @export
params_dialogue_hlm <- function(
  min_cells_per_sample = 2L,
  use_tme_qc = TRUE,
  use_cell_quality = TRUE,
  satterthwaite = TRUE
) {
  # checks
  checkmate::qassert(min_cells_per_sample, "I1[0,)")
  checkmate::qassert(use_tme_qc, "B1")
  checkmate::qassert(use_cell_quality, "B1")
  checkmate::qassert(satterthwaite, "B1")

  list(
    min_cells_per_sample = min_cells_per_sample,
    use_tme_qc = use_tme_qc,
    use_cell_quality = use_cell_quality,
    satterthwaite = satterthwaite
  )
}

#' Wrapper function for the DIALOGUE refinement parameters
#'
#' @description
#' Stage three of DIALOGUE: the cross-partner meta-analysis that decides which
#' genes survive, and the non-negative refit of the programme scores onto them.
#'
#' @details
#' Two gene lists come out. The permissive one asks only for a Fisher-combined
#' p below `permissive_p`; the strict one is looser on the p-value but also
#' demands that *every* partner supports the gene. They are not nested by
#' threshold, they are nested by evidence, and the strict list is the one to
#' quote.
#'
#' @param support_p Numeric. Adjusted p below which one partner counts as
#' supporting a gene. Must be in `(0, 1]`.
#' @param min_support_fraction Numeric. Minimum supporting fraction for a
#' stratum to enter the staged fit. Must be in `[0, 1]`.
#' @param min_stratum Integer. Minimum genes in a stratum before it is worth
#' fitting.
#' @param early_stop_cor Numeric. Correlation between the original score and
#' the running fit at which the staged fit stops early. Must be in `(0, 1]`.
#' @param permissive_p Numeric. Fisher-combined p for the permissive gene list,
#' where a gene is carried by partner support rather than by a positive
#' coefficient. Must be in `(0, 1]`.
#' @param strict_p Numeric. Fisher-combined p for the strict gene list, which
#' also demands that every partner supports the gene. Must be in `(0, 1]`.
#'
#' @returns A list with the stage three DIALOGUE parameters.
#'
#' @references Jerby-Arnon & Regev, Nature Biotechnology, 2022
#'
#' @export
params_dialogue_refine <- function(
  support_p = 0.1,
  min_support_fraction = 1 / 3,
  min_stratum = 5L,
  early_stop_cor = 0.95,
  permissive_p = 1e-3,
  strict_p = 0.05
) {
  # checks
  checkmate::qassert(support_p, "N1(0,1]")
  checkmate::qassert(min_support_fraction, "N1[0,1]")
  checkmate::qassert(min_stratum, "I1[0,)")
  checkmate::qassert(early_stop_cor, "N1(0,1]")
  checkmate::qassert(permissive_p, "N1(0,1]")
  checkmate::qassert(strict_p, "N1(0,1]")

  list(
    support_p = support_p,
    min_support_fraction = min_support_fraction,
    min_stratum = min_stratum,
    early_stop_cor = early_stop_cor,
    permissive_p = permissive_p,
    strict_p = strict_p
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
#' `c("nndescent", "ivf", "hnsw", "annoy", "exhaustive")`. Defaults to
#' `"nndescent"`. `"exhaustive"` scans every cell, so it returns the true
#' nearest neighbour rather than an approximation, at a cost that grows with
#' the number of cells.
#' @param knn List. Optional overrides for kNN parameters. See
#' [bixverse::params_knn_defaults()] for available parameters: `k`,
#' `knn_method`, `ann_dist`, `search_budget`, `n_trees`, `delta`,
#' `diversify_prob`, `ef_budget`, `m`, `ef_construction`, `ef_search`, `n_list`
#' and `n_probe`.
#'
#' @returns A list with the MiloR parameters.
#'
#' @export
params_sc_miloR <- function(
  prop = 0.2,
  k_refine = 20L,
  refinement_strategy = c("index", "approximate", "bruteforce"),
  index_type = c("nndescent", "ivf", "hnsw", "annoy", "exhaustive"),
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
#' updates. Defaults to `TRUE`.
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
  pruning = TRUE,
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
#' @param k_ith Optional integer. The k-ith neighbour to use for the kernel.
#' Defaults to `k %/% 2`.
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

  res <- c(
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
  # for easier detection down the line
  class(res) <- c("params_sc_harmony", "list")
  res
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

  res <- c(
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

  # for some tricks with symphony
  class(res) <- c("params_sc_harmony_v2", "list")

  res
}

### seurat CCA -----------------------------------------------------------------

#' Wrapper function for the Seurat CCA parameters
#'
#' @param num_cc Integer. Number of canonical correlation dimensions to compute
#' for the anchor space. Defaults to `30L`. The effective rank used is
#' `max(num_cc, dims)`.
#' @param dims Integer. Number of dimensions used for the anchor kNN queries
#' and the size of the returned embedding. Defaults to `30L`.
#' @param k_anchor Integer. Neighbourhood size for the mutual nearest neighbour
#' anchor search. Defaults to `5L`.
#' @param k_filter Integer. Neighbourhood size for the gene-space anchor filter.
#' Defaults to `200L`.
#' @param k_score Integer. Neighbourhood size for the shared-neighbour anchor
#' scoring. Defaults to `30L`.
#' @param k_weight Integer. Neighbourhood size for the kernel weights applied
#' during the correction. Defaults to `100L`.
#' @param n_top_features Integer. Number of top-loading genes used for the
#' gene-space anchor filter. Defaults to `200L`.
#' @param l2_norm Boolean. Shall the canonical correlation embedding be
#' L2-normalised per cell. Defaults to `TRUE`.
#' @param sd Numeric. Bandwidth divisor of the Gaussian kernel used for the
#' anchor weights. Defaults to `1.0`.
#' @param knn List. Optional overrides for kNN parameters. See
#' [bixverse::params_knn_defaults()] for available parameters: `k`,
#' `knn_method`, `ann_dist`, `search_budget`, `n_trees`, `delta`,
#' `diversify_prob`, `ef_budget`, `m`, `ef_construction`, `ef_search`, `n_list`
#' and `n_probe`. Note that `k` is unused here, the neighbourhood sizes come
#' from `k_anchor`, `k_filter`, `k_score` and `k_weight`.
#' @param pca Named list. Parameters to feed through to the optional
#' recalculation of the PCA, see [params_sc_pca()].
#'
#' @returns A list with the Seurat CCA parameters.
#'
#' @export
#'
#' @references Stuart, et al., Cell, 2019
params_sc_seurat_cca <- function(
  num_cc = 30L,
  dims = 30L,
  k_anchor = 5L,
  k_filter = 200L,
  k_score = 30L,
  k_weight = 100L,
  n_top_features = 200L,
  l2_norm = TRUE,
  sd = 1.0,
  knn = list(),
  pca = params_sc_pca()
) {
  checkmate::qassert(num_cc, "I1[1,)")
  checkmate::qassert(dims, "I1[1,)")
  checkmate::qassert(k_anchor, "I1[1,)")
  checkmate::qassert(k_filter, "I1[0,)")
  checkmate::qassert(k_score, "I1[1,)")
  checkmate::qassert(k_weight, "I1[1,)")
  checkmate::qassert(n_top_features, "I1[1,)")
  checkmate::qassert(l2_norm, "B1")
  checkmate::qassert(sd, "N1(0,)")

  knn_params <- modifyList(
    params_knn_defaults(),
    modifyList(list(ann_dist = "cosine"), knn, keep.null = TRUE),
    keep.null = TRUE
  )

  c(
    list(
      num_cc = num_cc,
      dims = dims,
      k_anchor = k_anchor,
      k_filter = k_filter,
      k_score = k_score,
      k_weight = k_weight,
      n_top_features = n_top_features,
      l2_norm = l2_norm,
      sd = sd
    ),
    knn_params,
    pca
  )
}

### seurat rPCA ----------------------------------------------------------------

#' Wrapper function for the Seurat rPCA parameters
#'
#' @param dims Integer. Number of dimensions used for the per-batch PCA
#' projections, the anchor kNN queries and the size of the returned embedding.
#' Defaults to `30L`.
#' @param k_anchor Integer. Neighbourhood size for the mutual nearest neighbour
#' anchor search. Defaults to `5L`.
#' @param k_score Integer. Neighbourhood size for the shared-neighbour anchor
#' scoring. Defaults to `30L`.
#' @param k_weight Integer. Neighbourhood size for the kernel weights applied
#' during the correction. Defaults to `100L`.
#' @param l2_norm Boolean. Shall the projected embeddings be L2-normalised per
#' cell. Defaults to `TRUE`.
#' @param sd Numeric. Bandwidth divisor of the Gaussian kernel used for the
#' anchor weights. Defaults to `1.0`.
#' @param knn List. Optional overrides for kNN parameters. See
#' [bixverse::params_knn_defaults()] for available parameters: `k`,
#' `knn_method`, `ann_dist`, `search_budget`, `n_trees`, `delta`,
#' `diversify_prob`, `ef_budget`, `m`, `ef_construction`, `ef_search`, `n_list`
#' and `n_probe`. Note that `k` is unused here, the neighbourhood sizes come
#' from `k_anchor`, `k_score` and `k_weight`.
#' @param pca Named list. Parameters to feed through to the optional
#' recalculation of the PCA, see [params_sc_pca()].
#'
#' @returns A list with the Seurat rPCA parameters. Note that rPCA has no
#' `num_cc`, `k_filter` or `n_top_features`, the gene-space anchor filter is
#' CCA-only in Seurat.
#'
#' @export
#'
#' @references Stuart, et al., Cell, 2019
params_sc_seurat_rpca <- function(
  dims = 30L,
  k_anchor = 5L,
  k_score = 30L,
  k_weight = 100L,
  l2_norm = TRUE,
  sd = 1.0,
  knn = list(),
  pca = params_sc_pca()
) {
  checkmate::qassert(dims, "I1[1,)")
  checkmate::qassert(k_anchor, "I1[1,)")
  checkmate::qassert(k_score, "I1[1,)")
  checkmate::qassert(k_weight, "I1[1,)")
  checkmate::qassert(l2_norm, "B1")
  checkmate::qassert(sd, "N1(0,)")

  knn_params <- modifyList(
    params_knn_defaults(),
    modifyList(list(ann_dist = "cosine"), knn, keep.null = TRUE),
    keep.null = TRUE
  )

  c(
    list(
      dims = dims,
      k_anchor = k_anchor,
      k_score = k_score,
      k_weight = k_weight,
      l2_norm = l2_norm,
      sd = sd
    ),
    knn_params,
    pca
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

## trajectory ------------------------------------------------------------------

### palantir -------------------------------------------------------------------

#' Wrapper function for Palantir parameters
#'
#' @description
#' Parameters controlling the Palantir trajectory inference. The kNN graph you
#' hand to [bixverse::run_palantir_sc()] feeds the diffusion kernel. The
#' geodesics are measured over a second kNN graph that Palantir builds
#' internally on the multiscale space, and `knn` is what controls that one. It
#' is a different knob from `k` in the kNN parameter block, which sizes the
#' backend index. Palantir overrides `k` and `ann_dist` for its internal search,
#' so only `knn_method` and the backend tuning parameters have an effect.
#'
#' @param n_dcs Integer. Diffusion components to extract before the multiscale
#' scaling. Defaults to `10L`.
#' @param n_eigs Optional integer. Eigenvectors to retain, not components: the
#' trivial leading eigenvector is counted here and then dropped, so `3L` leaves
#' two multiscale components. If `NULL`, the count is picked from the largest
#' eigengap, as the reference does.
#' @param knn Integer. Neighbours for the geodesic graph over the multiscale
#' space, in the reference's self-inclusive convention. Defaults to `30L`.
#' @param num_waypoints Integer. Target waypoint count for the max-min sampler.
#' Defaults to `1200L`.
#' @param scale_components Boolean. Min-max scale each multiscale component to
#' `[0, 1]` before any distance is taken. Defaults to `TRUE`.
#' @param use_early_cell_as_start Boolean. Use the provided early cell directly
#' rather than snapping it to the nearest diffusion-map boundary cell. Defaults
#' to `TRUE`, which deviates from the reference: the boundary candidate set is
#' at most two cells per multiscale component, so on a branching manifold the
#' snap can move a root cell onto a branch tip and run the trajectory backwards.
#' @param max_iterations Integer. Iteration cap for the pseudotime refinement.
#' Defaults to `25L`.
#' @param branch_prob_threshold Numeric. Fate probabilities below this are
#' zeroed. Defaults to `0.01`.
#' @param lanczos_basis_size Optional integer. Krylov basis vectors held at once
#' during the diffusion-map eigendecomposition. If `NULL`, derived from the
#' requested component count.
#' @param lanczos_max_restarts Integer. Maximum restart cycles for the Lanczos
#' solver. Defaults to `16L`.
#' @param lanczos_tol Numeric. Relative residual tolerance for the Lanczos
#' solver. Defaults to `1e-8`.
#' @param knn_params List. Optional overrides for the kNN parameters of the
#' internal multiscale search. See [bixverse::params_knn_defaults()] for
#' available parameters: `k`, `knn_method`, `ann_dist`, `search_budget`,
#' `n_trees`, `delta`, `diversify_prob`, `ef_budget`, `m`, `ef_construction`,
#' `ef_search`, `n_list` and `n_probe`.
#'
#' @returns A named flat list with all Palantir parameters.
#'
#' @export
#'
#' @references Setty, et al., Nat. Biotechnol., 2019.
params_sc_palantir <- function(
  n_dcs = 10L,
  n_eigs = NULL,
  knn = 30L,
  num_waypoints = 1200L,
  scale_components = TRUE,
  use_early_cell_as_start = TRUE,
  max_iterations = 25L,
  branch_prob_threshold = 0.01,
  lanczos_basis_size = NULL,
  lanczos_max_restarts = 16L,
  lanczos_tol = 1e-8,
  knn_params = list()
) {
  # checks
  checkmate::qassert(n_dcs, "I1[3,)")
  checkmate::qassert(n_eigs, c("0", "I1[3,)"))
  checkmate::qassert(knn, "I1[6,)")
  checkmate::qassert(num_waypoints, "I1[1,)")
  checkmate::qassert(scale_components, "B1")
  checkmate::qassert(use_early_cell_as_start, "B1")
  checkmate::qassert(max_iterations, "I1[2,)")
  checkmate::qassert(branch_prob_threshold, "N1[0,1]")
  checkmate::qassert(lanczos_basis_size, c("0", "I1[1,)"))
  checkmate::qassert(lanczos_max_restarts, "I1[1,)")
  checkmate::qassert(lanczos_tol, "N1(0,)")

  knn_defaults <- modifyList(
    params_knn_defaults(),
    knn_params,
    keep.null = TRUE
  )

  c(
    list(
      n_dcs = n_dcs,
      n_eigs = n_eigs,
      knn = knn,
      num_waypoints = num_waypoints,
      scale_components = scale_components,
      use_early_cell_as_start = use_early_cell_as_start,
      max_iterations = max_iterations,
      branch_prob_threshold = branch_prob_threshold,
      lanczos_basis_size = lanczos_basis_size,
      lanczos_max_restarts = lanczos_max_restarts,
      lanczos_tol = lanczos_tol
    ),
    knn_defaults
  )
}

### magic ----------------------------------------------------------------------

#' Wrapper function for MAGIC imputation parameters
#'
#' @description
#' Parameters controlling the MAGIC imputation run by
#' [bixverse::run_magic_sc()]. The defaults mirror the reference
#' implementation.
#'
#' @param n_steps Integer. Diffusion steps applied to the counts. Defaults to
#' `3L`. Zero is legal and hands back the un-imputed values, which is a cheap
#' way to compare the two.
#' @param clip_threshold Numeric. Imputed values below this are zeroed after
#' the last step. Defaults to `0.01`.
#' @param gene_batch_size Integer. Genes streamed off the binary store per
#' block. Bounds the scratch memory and is clamped to the number of requested
#' genes. Defaults to `1000L`.
#' @param layer String. One of `c("norm", "raw")`. Which stored layer to
#' impute. The operator preserves per-cell mass, so imputed values sit on the
#' scale of whatever went in: imputing raw counts and imputing log-normalised
#' counts are different operations rather than the same one rescaled. Defaults
#' to `"norm"`.
#' @param allow_large Boolean. Skip the output size guard. The dense output is
#' capped at 1e9 elements, i.e. 4 GB of `f32`. Defaults to `FALSE`.
#'
#' @returns A named flat list with all MAGIC parameters.
#'
#' @export
#'
#' @references van Dijk, et al., Cell, 2018.
params_sc_magic <- function(
  n_steps = 3L,
  clip_threshold = 0.01,
  gene_batch_size = 1000L,
  layer = c("norm", "raw"),
  allow_large = FALSE
) {
  layer <- match.arg(layer)

  # checks
  checkmate::qassert(n_steps, "I1[0,)")
  checkmate::qassert(clip_threshold, "N1[0,)")
  checkmate::qassert(gene_batch_size, "I1[1,)")
  checkmate::assertChoice(layer, c("norm", "raw"))
  checkmate::qassert(allow_large, "B1")

  list(
    n_steps = n_steps,
    clip_threshold = clip_threshold,
    gene_batch_size = gene_batch_size,
    layer = layer,
    allow_large = allow_large
  )
}

### gene trends ----------------------------------------------------------------

#' Wrapper function for the branch cell selection parameters
#'
#' @description
#' Parameters controlling which cells [bixverse::run_gene_trends_sc()] assigns
#' to each branch. The threshold on a fate's probability is an expanding
#' quantile over the pseudotime-sorted cells, made monotone with a cumulative
#' maximum, so a fate's bar can only rise as differentiation proceeds. The
#' defaults are the reference ones.
#'
#' @param q Numeric. Upper-tail quantile of the fate probability used as the
#' threshold. Defaults to `0.01`.
#' @param eps Numeric. Slack subtracted from the threshold before the
#' comparison. Defaults to `0.01`.
#' @param resolution Integer. Number of pseudotime buckets, capped at the cell
#' count. Defaults to `500L`.
#'
#' @returns A named flat list with all branch selection parameters.
#'
#' @export
#'
#' @references Setty, et al., Nat. Biotechnol., 2019.
params_sc_branch_selection <- function(
  q = 0.01,
  eps = 0.01,
  resolution = 500L
) {
  # checks
  checkmate::qassert(q, "N1[0,1]")
  checkmate::qassert(eps, "N1[0,)")
  checkmate::qassert(resolution, "I1[1,)")

  list(
    q = q,
    eps = eps,
    resolution = resolution
  )
}

#' Wrapper function for gene trend parameters
#'
#' @description
#' Parameters controlling the landmark Gaussian process that
#' [bixverse::run_gene_trends_sc()] fits per branch. The kernel is a
#' Matern-5/2 one and the prediction grid doubles as the landmark set.
#'
#' The defaults come from the reference and are prior-dominated. Palantir's
#' pseudotime is min-max scaled to `[0, 1]`, so a `length_scale` of `1.0` spans
#' the entire domain and a `sigma` of `1.0` sits at roughly the signal scale of
#' log-normalised expression. The posterior will flatten genuine transient
#' structure and resolve almost any gene into a smooth monotone or
#' single-peaked curve. That is a presentation choice, not inference. Shorten
#' `length_scale` before believing a bump.
#'
#' @param resolution Integer. Grid points per branch. Kept at the default even
#' when a branch holds fewer cells, as the reference does. Defaults to `500L`.
#' @param weighting String. One of `c("hard_mask", "fate_probability")`. With
#' `"hard_mask"` every selected cell enters its branch's fit with equal weight,
#' which is what the reference does. With `"fate_probability"` every cell
#' enters every fit weighted by its fate probability, which is more defensible:
#' a cell at 0.6 is not a member. Defaults to `"hard_mask"`.
#' @param length_scale Numeric. Matern-5/2 length scale. Defaults to `1.0`.
#' @param sigma Numeric. Noise standard deviation. Defaults to `1.0`.
#' @param jitter Numeric. Added to the landmark covariance diagonal before the
#' Cholesky. Defaults to `1e-6`.
#' @param max_jitter_retries Integer. Times the jitter is raised and the
#' Cholesky retried before giving up. Defaults to `3L`.
#' @param chunk_size Integer. Training points held at once when accumulating
#' the cross-covariance. Defaults to `2048L`.
#'
#' @returns A named flat list with all gene trend parameters. The Gaussian
#' process hyperparameters sit at the same level as the trend ones, not in a
#' nested block.
#'
#' @export
#'
#' @references Setty, et al., Nat. Biotechnol., 2019.
params_sc_gene_trends <- function(
  resolution = 500L,
  weighting = c("hard_mask", "fate_probability"),
  length_scale = 1.0,
  sigma = 1.0,
  jitter = 1e-6,
  max_jitter_retries = 3L,
  chunk_size = 2048L
) {
  weighting <- match.arg(weighting)

  # checks
  checkmate::qassert(resolution, "I1[2,)")
  checkmate::assertChoice(weighting, c("hard_mask", "fate_probability"))
  checkmate::qassert(length_scale, "N1(0,)")
  checkmate::qassert(sigma, "N1(0,)")
  checkmate::qassert(jitter, "N1[0,)")
  checkmate::qassert(max_jitter_retries, "I1[0,)")
  checkmate::qassert(chunk_size, "I1[1,)")

  list(
    resolution = resolution,
    weighting = weighting,
    length_scale = length_scale,
    sigma = sigma,
    jitter = jitter,
    max_jitter_retries = max_jitter_retries,
    chunk_size = chunk_size
  )
}

## lda -------------------------------------------------------------------------

#' Wrapper function for the LDA parameters
#'
#' @description
#' Solver options for the variational Bayes latent Dirichlet allocation, see
#' [bixverse::run_lda()].
#'
#' @details
#' The defaults follow pycisTopic, so the knobs mean the same thing on both
#' sides. `alpha_by_topic = TRUE` turns `alpha` into the Griffiths and
#' Steyvers `50 / k` heuristic that cisTopic defaults to; set it to `FALSE` if
#' you want `alpha` taken literally.
#'
#' `learning = "batch"` sweeps every document once per iteration and is
#' monotone in the bound. `"online"` takes decaying steps from shuffled
#' mini-batches, which reaches a usable fit in far fewer passes on a large
#' corpus at the cost of that guarantee. `batch_size` and `n_epochs` are only
#' read by the online variant.
#'
#' @param alpha Float. Dirichlet prior on the document-topic distributions.
#' @param alpha_by_topic Boolean. Shall `alpha` be divided by the topic count.
#' @param eta Float. Dirichlet prior on the topic-term distributions.
#' @param eta_by_topic Boolean. Shall `eta` be divided by the topic count.
#' @param max_iter Integer. Maximum outer iterations. Ignored by the online
#' variant, which counts epochs instead.
#' @param tol Float. Relative change in the bound below which the solver stops.
#' @param inner_max_iter Integer. Maximum fixed-point iterations of the
#' per-document E-step.
#' @param inner_tol Float. Relative L1 change in the variational parameters
#' below which the per-document E-step stops.
#' @param check_every Integer. Iterations between bound evaluations.
#' @param learning String. One of `c("batch", "online")`.
#' @param batch_size Integer. Documents per mini-batch. Online only.
#' @param n_epochs Integer. Passes over the corpus. Online only.
#'
#' @returns A list with the LDA parameters.
#'
#' @references Hoffman, Blei and Bach, NIPS, 2010; Bravo Gonzalez-Blas, et al.,
#' Nat Methods, 2019
#'
#' @export
params_lda <- function(
  alpha = 50.0,
  alpha_by_topic = TRUE,
  eta = 0.1,
  eta_by_topic = FALSE,
  max_iter = 150L,
  tol = 1e-3,
  inner_max_iter = 100L,
  inner_tol = 1e-3,
  check_every = 10L,
  learning = c("batch", "online"),
  batch_size = 1024L,
  n_epochs = 10L
) {
  learning <- match.arg(learning)

  # checks
  checkmate::qassert(alpha, "N1(0,)")
  checkmate::qassert(alpha_by_topic, "B1")
  checkmate::qassert(eta, "N1(0,)")
  checkmate::qassert(eta_by_topic, "B1")
  checkmate::qassert(max_iter, "I1[1,)")
  checkmate::qassert(tol, "N1(0,)")
  checkmate::qassert(inner_max_iter, "I1[1,)")
  checkmate::qassert(inner_tol, "N1(0,)")
  checkmate::qassert(check_every, "I1[1,)")
  checkmate::assertChoice(learning, c("batch", "online"))
  checkmate::qassert(batch_size, "I1[1,)")
  checkmate::qassert(n_epochs, "I1[1,)")

  list(
    alpha = alpha,
    alpha_by_topic = alpha_by_topic,
    eta = eta,
    eta_by_topic = eta_by_topic,
    max_iter = max_iter,
    tol = tol,
    inner_max_iter = inner_max_iter,
    inner_tol = inner_tol,
    check_every = check_every,
    learning = learning,
    batch_size = batch_size,
    n_epochs = n_epochs
  )
}

## differential expression -----------------------------------------------------

### edgeR quasi-likelihood -----------------------------------------------------

#' Wrapper function for parameters for the edgeR quasi-likelihood workflow
#'
#' @description
#' Parameters for the edgeR quasi-likelihood chain, implemented in Rust via the
#' `edge-rs` crate and gated against edgeR 4.8.2. Defaults are edgeR's own.
#'
#' The `legacy` switch picks between two genuinely different pipelines. The
#' current route estimates its own dispersion from the most abundant genes and
#' skips `estimateDisp()`, which is where most of the runtime went and is
#' edgeR 4's own recommendation. The legacy route shrinks the raw residual
#' deviance, needs a dispersion handed to it, and is the only one where the
#' Poisson bound bites.
#'
#' @param norm_method String. Library size normalisation. One of
#' `c("TMM", "TMMwsp", "RLE", "upperquartile", "none")`. Defaults to `"TMM"`.
#' `"none"` leaves every factor at one, which is what Milo's `logMS` amounts
#' to.
#' @param filter Boolean. Run `filterByExpr()` before fitting. Defaults to
#' `TRUE`. Turn this off for anything that is not gene expression, e.g. Milo
#' neighbourhood counts, where the heuristic means nothing.
#' @param min_mean Numeric. Drop features whose mean count across samples is
#' below this. Applied on top of `filter`. Defaults to `0` (off).
#' @param robust Boolean. Robust empirical Bayes squeezing, giving outlier
#' features their own smaller prior degrees of freedom. Defaults to `FALSE`.
#' @param legacy Boolean. Take edgeR's pre-4.0 quasi-likelihood pipeline.
#' Defaults to `FALSE`.
#'
#' @returns A list with the edgeR quasi-likelihood parameters.
#'
#' @references Chen, Lun and Smyth, F1000Research, 2016
#'
#' @export
params_edger_ql <- function(
  norm_method = c("TMM", "TMMwsp", "RLE", "upperquartile", "none"),
  filter = TRUE,
  min_mean = 0,
  robust = FALSE,
  legacy = FALSE
) {
  norm_method <- match.arg(norm_method)

  # checks
  checkmate::assertChoice(
    norm_method,
    c("TMM", "TMMwsp", "RLE", "upperquartile", "none")
  )
  checkmate::qassert(filter, "B1")
  checkmate::qassert(min_mean, "N1[0,)")
  checkmate::qassert(robust, "B1")
  checkmate::qassert(legacy, "B1")

  list(
    norm_method = norm_method,
    filter = filter,
    min_mean = min_mean,
    robust = robust,
    legacy = legacy
  )
}
