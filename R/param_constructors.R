# param constructors -----------------------------------------------------------

## defaults --------------------------------------------------------------------

#' Helper function to generate kNN defaults
#'
#' @returns A list with default parameters for kNN searches. Following
#' parameters:
#' \itemize{
#'  \item k - Number of neighbours. Defaults to `15L`.
#'  \item knn_method - Which of method to use for the approximate nearest
#'  neighbour search. Defaults to `"annoy"`. The implementations are:
#'  `c("annoy", "hnsw", "nndescent")`.
#'  \item ann_dist - Which distance metric to use for the approximate nearest
#'  neighbour search. Defaults to `"euclidean"`. The implementations are
#'  `c("euclidean", "cosine")`.
#'  \item search_budget - Search budget per tree for Annoy. Defaults to `100L`.
#'  \item n_trees - Number of trees to generate for Annoy. Defaults to `100L`.
#'  \item nn_max_iter - Maximum iterations for NNDescent. Defaults to `15L`.
#'  \item rho - Sampling rate for NNDescent. Defaults to `1.0`.
#'  \item delta - Early termination criterium for NNDescent. Defaults to
#'  `0.001`.
#' }
#'
#' @export
params_knn_defaults <- function() {
  list(
    k = 15L,
    knn_method = "annoy",
    ann_dist = "euclidean",
    search_budget = 100L,
    n_trees = 100L,
    nn_max_iter = 15L,
    rho = 1.0,
    delta = 0.001
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
#' `knn_method`, `dist_metric`, `search_budget`, `n_trees`, `nn_max_iter`,
#' `rho`, `delta`.
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
  knn = list(k = 0L, knn_method = "hnsw")
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
#' `knn_method`, `dist_metric`, `search_budget`, `n_trees`, `nn_max_iter`,
#' `rho`, `delta`.
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
  knn = list(k = 0L)
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

## neighbours ------------------------------------------------------------------

#' Wrapper function for parameters for neighbour identification in single cell
#'
#' @param full_snn Boolean. Shall the full shared nearest neighbour graph
#' be generated that generates edges between all cells instead of between
#' only neighbours.
#' @param pruning Numeric. Weights below this threshold will be set to 0 in
#' the generation of the sNN graph.
#' @param snn_similarity String. One of `c("rank", "jaccard")`. The Jaccard
#' similarity calculates the Jaccard between the neighbours, whereas the rank
#' method calculates edge weights based on the ranking of shared neighbours.
#' For the rank method, the weight is determined by finding the shared
#' neighbour with the lowest combined rank across both cells, where
#' lower-ranked (closer) shared neighbours result in higher edge weights
#' Both methods produce weights normalised to the range `[0, 1]`.
#' @param knn List. Optional overrides for kNN parameters. See
#' [bixverse::params_knn_defaults()] for available parameters: `k`,
#' `knn_method`, `ann_dist`, `search_budget`, `n_trees`, `nn_max_iter`,
#' `rho`, `delta`. Note: This function uses `k = 15L`, `ann_dist = "cosine"`,
#' and `nn_max_iter = 25L` as defaults.
#'
#' @returns A list with the neighbour parameters.
#'
#' @export
params_sc_neighbours <- function(
  full_snn = FALSE,
  pruning = 1 / 15,
  snn_similarity = c("rank", "jaccard"),
  knn = list()
) {
  snn_similarity <- match.arg(snn_similarity)

  checkmate::qassert(full_snn, "B1")
  checkmate::qassert(pruning, "N1[0, 1]")
  checkmate::assertChoice(snn_similarity, c("rank", "jaccard"))

  knn_params <- modifyList(
    params_knn_defaults(),
    modifyList(
      list(k = 15L, ann_dist = "cosine"),
      knn,
      keep.null = TRUE
    ),
    keep.null = TRUE
  )

  checkmate::qassert(knn_params$k, "I1")
  checkmate::assertChoice(
    knn_params$knn_method,
    c("annoy", "hnsw", "nndescent")
  )
  checkmate::assertChoice(knn_params$ann_dist, c("cosine", "euclidean"))
  checkmate::qassert(knn_params$n_trees, "I1")
  checkmate::qassert(knn_params$search_budget, "N1")
  checkmate::qassert(knn_params$nn_max_iter, "I1")
  checkmate::qassert(knn_params$rho, "N1")
  checkmate::qassert(knn_params$delta, "N1")

  c(
    list(
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
#' `knn_method`, `ann_dist`, `search_budget`, `n_trees`, `nn_max_iter`,
#' `rho`, `delta`. Note: This function uses `k = 15L` and `nn_max_iter = 25L`
#' as defaults.
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

  checkmate::qassert(knn_params$k, "I1")
  checkmate::assertChoice(
    knn_params$knn_method,
    c("annoy", "hnsw", "nndescent")
  )
  checkmate::qassert(knn_params$n_trees, "I1")
  checkmate::qassert(knn_params$search_budget, "N1")
  checkmate::qassert(knn_params$nn_max_iter, "I1")
  checkmate::qassert(knn_params$rho, "N1")
  checkmate::qassert(knn_params$delta, "N1")

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
#' `c("danb", "bernoulli", "normal")`. Defaults to `"normal"`.
#' @param normalise Boolean. Shall the data be normalised. Defaults to `TRUE`.
#' @param knn List. Optional overrides for kNN parameters. See
#' [bixverse::params_knn_defaults()] for available parameters: `k`,
#' `knn_method`, `ann_dist`, `search_budget`, `n_trees`, `nn_max_iter`,
#' `rho`, `delta`. Note: Parameters are mapped to Rust structure names
#' (`n_trees` becomes `n_tree`, `nn_max_iter` becomes `max_iter`).
#'
#' @returns A list with the HotSpot parameters.
#'
#' @export
params_sc_hotspot <- function(
  model = c("normal", "danb", "bernoulli"),
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

  checkmate::qassert(knn_params$k, "I1")
  checkmate::assertChoice(
    knn_params$knn_method,
    c("annoy", "hnsw", "nndescent")
  )
  checkmate::assertChoice(knn_params$ann_dist, c("cosine", "euclidean"))
  checkmate::qassert(knn_params$n_trees, "I1")
  checkmate::qassert(knn_params$search_budget, "I1")
  checkmate::qassert(knn_params$nn_max_iter, "I1")
  checkmate::qassert(knn_params$rho, "N1")
  checkmate::qassert(knn_params$delta, "N1")

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

## metacells -------------------------------------------------------------------

### meta cell (hdWGCNA) --------------------------------------------------------

#' Wrapper function for parameters for meta cell generation
#'
#' @param max_shared Integer. Maximum number of allowed shared neighbours for
#' the meta cell to be considered. Defaults to `15L`.
#' @param target_no_metacells Integer. Target number of meta-cells to generate.
#' Defaults to `1000L`.
#' @param max_iter Integer. Maximum number of iterations for the algorithm.
#' Defaults to `5000L`.
#' @param knn List. Optional overrides for kNN parameters. See
#' [bixverse::params_knn_defaults()] for available parameters: `k`,
#' `knn_method`, `ann_dist`, `search_budget`, `n_trees`, `nn_max_iter`,
#' `rho`, `delta`. Note: This function uses `k = 25L` and `ann_dist = "cosine"`
#' as defaults.
#'
#' @returns A list with the metacell parameters.
#'
#' @export
params_sc_metacells <- function(
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

  checkmate::qassert(knn_params$k, "I1")
  checkmate::assertChoice(
    knn_params$knn_method,
    c("annoy", "hnsw", "nndescent")
  )
  checkmate::assertChoice(knn_params$ann_dist, c("cosine", "euclidean"))
  checkmate::qassert(knn_params$n_trees, "I1")
  checkmate::qassert(knn_params$search_budget, "I1")
  checkmate::qassert(knn_params$nn_max_iter, "I1")
  checkmate::qassert(knn_params$rho, "N1")
  checkmate::qassert(knn_params$delta, "N1")

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
#' @param knn List. Optional overrides for kNN parameters. See
#' [bixverse::params_knn_defaults()] for available parameters: `k`,
#' `knn_method`, `ann_dist`, `search_budget`, `n_trees`, `nn_max_iter`,
#' `rho`, `delta`. Note: This function uses `k = 25L` as default.
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

  knn_params <- modifyList(
    params_knn_defaults(),
    modifyList(list(k = 25L), knn, keep.null = TRUE),
    keep.null = TRUE
  )

  checkmate::qassert(knn_params$k, "I1")
  checkmate::assertChoice(
    knn_params$knn_method,
    c("annoy", "hnsw", "nndescent")
  )
  checkmate::qassert(knn_params$n_trees, "I1")
  checkmate::qassert(knn_params$search_budget, "I1")
  checkmate::qassert(knn_params$nn_max_iter, "I1")
  checkmate::qassert(knn_params$rho, "N1")
  checkmate::qassert(knn_params$delta, "N1")

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
      pruning_threshold = pruning_threshold
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
#' @param linkage_dist String. Which type of distance metric to use for the
#' linkage. Defaults to `"average"`.
#' @param knn List. Optional overrides for kNN parameters. See
#' [bixverse::params_knn_defaults()] for available parameters: `k`,
#' `knn_method`, `ann_dist`, `search_budget`, `n_trees`, `nn_max_iter`,
#' `rho`, `delta`. Note: This function uses `k = 5L` and `ann_dist = "cosine"`
#' as defaults.
#'
#' @returns A list with the SuperCell parameters.
#'
#' @export
params_sc_supercell <- function(
  walk_length = 3L,
  graining_factor = 20.0,
  linkage_dist = c("complete", "average"),
  knn = list()
) {
  linkage_dist <- match.arg(linkage_dist)

  checkmate::qassert(walk_length, "I1")
  checkmate::qassert(graining_factor, "N1")
  checkmate::assertChoice(linkage_dist, c("complete", "average"))

  knn_params <- modifyList(
    params_knn_defaults(),
    modifyList(list(k = 5L, ann_dist = "cosine"), knn, keep.null = TRUE),
    keep.null = TRUE
  )

  checkmate::qassert(knn_params$k, "I1")
  checkmate::assertChoice(
    knn_params$knn_method,
    c("annoy", "hnsw", "nndescent")
  )
  checkmate::assertChoice(knn_params$ann_dist, c("cosine", "euclidean"))
  checkmate::qassert(knn_params$n_trees, "I1")
  checkmate::qassert(knn_params$search_budget, "I1")
  checkmate::qassert(knn_params$nn_max_iter, "I1")
  checkmate::qassert(knn_params$rho, "N1")
  checkmate::qassert(knn_params$delta, "N1")

  c(
    list(
      walk_length = walk_length,
      graining_factor = graining_factor,
      linkage_dist = linkage_dist
    ),
    knn_params
  )
}
