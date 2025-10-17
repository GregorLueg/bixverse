# param constructors -----------------------------------------------------------

## defaults --------------------------------------------------------------------

#' Helper function to generate kNN defaults
#'
#' @returns A list with default parameters for kNN searches. Following
#' parameters:
#' \itemize{
#'  \item k - Number of neighbours. Defaults to `15L`.
#'  \item knn_method - Which of method to use for the approximate nearest
#'  neighbour search. Defaults to `"annoy"`.
#'  \item dist_metric - Which distance metric to use for the approximate nearest
#'  neighbour search. Defaults to `"euclidean"`.
#'  \item search_budget - Search budget per tree for Annoy. Defaults to `100L`.
#'  \item n_trees - Number of trees to generate for Annoy. Defaults to `100L`.
#' }
#'
#' @export
params_knn_defaults <- function() {
  list(
    k = 0L,
    knn_method = "annoy",
    dist_metric = "euclidean",
    search_budget = 100L,
    n_trees = 100L
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
#'  \item loess_span - In case of
#'  \item search_budget - Search budget per tree for Annoy. Defaults to `100L`.
#'  \item n_trees - Number of trees to generate for Annoy. Defaults to `100L`.
#' }
#'
#' @export
params_hvg_defaults <- function() {
  list(
    min_gene_var_pctl = 0.7,
    hvg_method = "vst",
    loess_span = 0.3,
    clip_max = NULL
  )
}

#' Helper function to generate normalisation defaults for doublet detection.
#'
#' @return A list with the following parameters for normalisation specifically
#' designed for doublet detection methods:
#' \itemize{
#'  \item log_transform - Boolean. Shall the counts be log-normalised.
#'  Defaults to `TRUE`.
#'  \item mean_center - Boolean. Shall mean centreing be applied. Defaults
#'  to `FALSE`.
#'  \item normalise_variance - Boolean. Shall the variance be normalised.
#'  Defaults to `FALSE`.
#'  \item target_size - Target library size. Defaults to `1e6`
#' }
#'
#' @export
params_norm_doublet_detection_defaults <- function() {
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
#' }
params_pca_defaults <- function() {
  list(
    no_pcs = 30L,
    random_svd = TRUE
  )
}

# constructors -----------------------------------------------------------------

## scrublet --------------------------------------------------------------------

#' Wrapper function for Scrublet doublet detection parameters
#'
#' @param sim_doublet_ratio Numeric. Number of doublets to simulate relative to
#' the number of observed cells. For example, 2.0 simulates twice as many
#' doublets as there are cells. Defaults to `1.5`.
#' @param expected_doublet_rate Numeric. Expected doublet rate for the
#' experiment, typically 0.05-0.10 depending on cell loading. Must be between
#' 0 and 1. Defaults to `0.1`.
#' @param stdev_doublet_rate Numeric. Uncertainty in the expected doublet rate.
#' Defaults to `0.02`.
#' @param n_bins Integer. Number of bins for histogram-based automatic threshold
#' detection. Typically 50-100. Defaults to `100L`.
#' @param manual_threshold Optional numeric. Manual doublet score threshold. If
#' `NULL` (default), threshold is automatically detected from simulated doublet
#' score distribution.
#' @param normalisation List. Optional overrides for normalisation parameters.
#' See [bixverse::params_norm_doublet_detection_defaults()] for available
#' parameters: `log_transform`, `mean_center`, `normalise_variance`,
#' `target_size`.
#' @param hvg List. Optional overrides for highly variable gene selection
#' parameters. See [bixverse::params_hvg_defaults()] for available parameters:
#' `min_gene_var_pctl`, `hvg_method`, `loess_span`, `clip_max`.
#' @param pca List. Optional overrides for PCA parameters. See
#' [bixverse::params_pca_defaults()] for available parameters: `no_pcs`,
#' `random_svd`.
#' @param knn List. Optional overrides for kNN parameters. See
#' [bixverse::params_knn_defaults()] for available parameters: `k`,
#' `knn_method`, `dist_metric`, `search_budget`, `n_trees`.
#'
#' @returns A named list with all Scrublet parameters, combining defaults with
#' any user-specified overrides.
#'
#' @export
params_scrublet <- function(
  sim_doublet_ratio = 1.5,
  expected_doublet_rate = 0.1,
  stdev_doublet_rate = 0.02,
  n_bins = 100L,
  manual_threshold = NULL,
  normalisation = list(),
  hvg = list(),
  pca = list(),
  knn = list()
) {
  # doublet simulation checks
  checkmate::qassert(sim_doublet_ratio, "N1(0,)")
  checkmate::qassert(expected_doublet_rate, "N1[0,1]")
  checkmate::qassert(stdev_doublet_rate, "N1[0,1]")
  checkmate::qassert(n_bins, "I1[10,)")
  if (!is.null(manual_threshold)) {
    checkmate::qassert(manual_threshold, "N1[0,)")
  }

  # generate final parameters
  params <- list(
    normalisation = modifyList(
      params_norm_doublet_detection_defaults(),
      normalisation,
      keep.null = TRUE
    ),
    hvg = modifyList(params_hvg_defaults(), hvg, keep.null = TRUE),
    pca = modifyList(params_pca_defaults(), pca, keep.null = TRUE),
    knn = modifyList(params_knn_defaults(), knn, keep.null = TRUE),
    sim_doublet_ratio = sim_doublet_ratio,
    expected_doublet_rate = expected_doublet_rate,
    stdev_doublet_rate = stdev_doublet_rate,
    n_bins = n_bins,
    manual_threshold = manual_threshold
  )

  params <- purrr::list_flatten(params, name_spec = "{inner}")

  params
}

## boost -----------------------------------------------------------------------

#' Wrapper function for Boost parameters
#'
#' @param boost_rate Numeric. Boosting rate for the algorithm. Must be between
#' 0 and 1. Defaults to `0.25`.
#' @param replace Boolean. Whether to use replacement during boosting. Defaults
#' to `FALSE`.
#' @param resolution Numeric. Resolution parameter for graph-based clustering.
#' Higher values lead to more clusters. Defaults to `1.0`.
#' @param n_iters Integer. Number of iterations to run the algorithm. Defaults
#' to `25L`.
#' @param p_thresh Numeric. P-value threshold for significance testing. Defaults
#' to `1e-7`.
#' @param voter_thresh Numeric. Voter threshold across iterations. Proportion of
#' iterations a cell must be assigned to a cluster to be considered a member.
#' Must be between 0 and 1. Defaults to `0.9`.
#' @param normalisation List. Optional overrides for normalisation parameters.
#' See [bixverse::params_norm_doublet_detection_defaults()] for available
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
#' `knn_method`, `dist_metric`, `search_budget`, `n_trees`.
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
  n_iters = 25L,
  p_thresh = 1e-7,
  voter_thresh = 0.9,
  normalisation = list(),
  hvg = list(),
  pca = list(),
  knn = list()
) {
  # checks
  checkmate::qassert(boost_rate, "N1[0,1]")
  checkmate::qassert(replace, "B1")
  checkmate::qassert(resolution, "N1(0,)")
  checkmate::qassert(n_iters, "I1[1,)")
  checkmate::qassert(p_thresh, "N1(0,)")
  checkmate::qassert(voter_thresh, "N1[0,1]")

  # generate final parameters
  params <- list(
    normalisation = modifyList(
      params_norm_doublet_detection_defaults(),
      normalisation,
      keep.null = TRUE
    ),
    hvg = modifyList(params_hvg_defaults(), hvg, keep.null = TRUE),
    pca = modifyList(params_pca_defaults(), pca, keep.null = TRUE),
    knn = modifyList(params_knn_defaults(), knn, keep.null = TRUE),
    boost_rate = boost_rate,
    replace = replace,
    resolution = resolution,
    n_iters = n_iters,
    p_thresh = p_thresh,
    voter_thresh = voter_thresh
  )

  params <- purrr::list_flatten(params, name_spec = "{inner}")

  assertScBoost(params)

  params
}
