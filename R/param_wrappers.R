# param wrapper ----------------------------------------------------------------

## ica -------------------------------------------------------------------------

#' Wrapper function for standard ICA parameters
#'
#' @param maxit Integer. Maximum number of iterations for ICA.
#' @param alpha Float. The alpha parameter for the logcosh version of ICA.
#' Should be between 1 to 2.
#' @param max_tol Numeric. Should be `0 < max_tol < 1`. Maximum tolerance of the
#' algorithm.
#' @param verbose Boolean. Controls verbosity of the function.
#'
#' @returns A list with the parameters for usage in subsequent functions.
#'
#' @export
params_ica_general <- function(
  maxit = 200L,
  alpha = 1.0,
  max_tol = 0.0001,
  verbose = FALSE
) {
  # Checks
  checkmate::qassert(maxit, "I1")
  checkmate::qassert(alpha, "R1[1, 2]")
  checkmate::qassert(max_tol, "R1(0, 1)")
  checkmate::qassert(verbose, "B1")
  # Params
  return(list(
    maxit = maxit,
    alpha = alpha,
    max_tol = max_tol,
    verbose = verbose
  ))
}


#' Wrapper function for ICA ncomp iterations
#'
#' @description Wrapper function to provide parameters through which ncomps to
#' iterate through.
#'
#' @param max_no_comp Integer. Maximum number of ncomp to test.
#' @param steps Integer. In which steps to move from 5 onwards.
#' @param custom_seq An integer vector. If you wish to provide a custom version
#' of no_comp to iterate through. If NULL, you will iterate through
#' `c(2, 3, 4, 5, 5 + step, ... max_no_comp - step, max_no_comp)`
#'
#' @returns A list with the parameters for usage in subsequent functions.
#'
#' @export
params_ica_ncomp <- function(max_no_comp = 75L, steps = 5L, custom_seq = NULL) {
  # Checks
  checkmate::qassert(max_no_comp, "I1")
  checkmate::qassert(steps, "I1")
  checkmate::qassert(custom_seq, c("0", "I+"))
  return(list(
    max_no_comp = max_no_comp,
    steps = steps,
    custom_seq = custom_seq
  ))
}

#' Wrapper function for ICA randomisation
#'
#'
#' @param cross_validate Boolean. Do you want to apply a cross-validation type
#' approach and split the data into `folds` folds to assess within data
#' stability of the component.
#' @param random_init Integer. Number of random initialisations to use.
#' @param folds Integer. Number of folds to use if `cross_validate` is set to
#' `TRUE`. To note, you will be running `random_init * folds` ICA runs.
#'
#' @returns A list with the parameters for usage in the subsequent functions.
#'
#' @export
params_ica_randomisation <- function(
  cross_validate = FALSE,
  random_init = 50L,
  folds = 10L
) {
  # Checks
  checkmate::qassert(cross_validate, "B1")
  checkmate::qassert(random_init, "I1")
  checkmate::qassert(folds, "I1")
  return(list(
    cross_validate = cross_validate,
    random_init = random_init,
    folds = folds
  ))
}

## graph-based stuff -----------------------------------------------------------

#' Wrapper function for graph generation
#'
#' @param epsilon Float. Defines the epsilon parameter for the radial basis
#' function. Defaults to 2, but should be ideally optimised.
#' @param min_cor Float. Minimum absolute correlation that needs to be
#' observed in either data set. Only relevant for differential correlation-based
#' graphs.
#' @param fdr_threshold Float. Maximum FDR for the differential correlation
#' p-value.
#' @param verbose Boolean. Controls verbosity of the graph generation function.
#'
#' @returns List with parameters for usage in subsequent function.
#'
#' @export
params_cor_graph <- function(
  epsilon = 2,
  min_cor = 0.2,
  fdr_threshold = 0.05,
  verbose = TRUE
) {
  # Checks
  checkmate::qassert(epsilon, "R1")
  checkmate::qassert(min_cor, "R1[0, 1]")
  checkmate::qassert(fdr_threshold, "R1[0, 1]")
  checkmate::qassert(verbose, "B1")
  # Return
  return(
    list(
      epsilon = epsilon,
      min_cor = min_cor,
      fdr_threshold = fdr_threshold,
      verbose = verbose
    )
  )
}

#' Wrapper function to generate resolution parameters for Leiden or Louvain
#' clustering.
#'
#' @param min_res Float. Minimum resolution to test.
#' @param max_res Float. Maximum resolution to test.
#' @param number_res Integer. Number of resolutions to test between the `max_res`
#' and `min_res.`
#'
#' @returns List with parameters for usage in subsequent function.
#'
#' @export
params_graph_resolution <- function(
  min_res = 0.1,
  max_res = 10,
  number_res = 15L
) {
  # Checks
  checkmate::qassert(min_res, "R1")
  checkmate::qassert(max_res, "R1")
  checkmate::qassert(number_res, "I1")
  # Return
  return(list(
    min_res = min_res,
    max_res = max_res,
    number_res = number_res
  ))
}

## community detections --------------------------------------------------------

#' Wrapper function to generate community detection parameters
#'
#' @param max_nodes Integer. Maximum number of nodes in a given community.
#' @param min_nodes Integer. Minimum number of nodes in a given community.
#' @param min_seed_nodes Integer. Minimum number of seed nodes within a
#' community.
#' @param initial_res Float. Initial resolution parameter to start with.
#' @param threshold_type String. One of `c("prop_based", "pval_based")`.
#' You can chose to include a certain proportion of the network with the highest
#' diffusion scores, or use p-values based on permutations.
#' @param network_threshold Float. The proportion of the network to include.
#' Used if `threshold_type = "prop_based"`.
#' @param pval_threshold Float. The maximum p-value for nodes to be included.
#' Used if `threshold_type = "pval_based"`.
#'
#' @returns List with parameters for usage in subsequent function.
#'
#' @export
params_community_detection <- function(
  max_nodes = 300L,
  min_nodes = 10L,
  min_seed_nodes = 2L,
  initial_res = 0.5,
  threshold_type = c("prop_based", "pval_based"),
  network_threshold = 0.5,
  pval_threshold = 0.1
) {
  threshold_type <- match.arg(threshold_type)

  # Checks
  checkmate::qassert(max_nodes, sprintf("I1[%i,)", min_nodes))
  checkmate::qassert(min_nodes, "I1")
  checkmate::qassert(min_seed_nodes, "I1")
  checkmate::qassert(initial_res, "N1")
  checkmate::assertChoice(threshold_type, c("prop_based", "pval_based"))
  checkmate::qassert(network_threshold, "N1(0, 1]")
  checkmate::qassert(pval_threshold, "N1(0, 1]")
  # Return
  return(
    list(
      max_nodes = max_nodes,
      min_nodes = min_nodes,
      min_seed_nodes = min_seed_nodes,
      initial_res = initial_res,
      threshold_type = threshold_type,
      network_threshold = network_threshold,
      pval_threshold = pval_threshold
    )
  )
}

## fgsea -----------------------------------------------------------------------

#' Wrapper function to generate GSEA parameters
#'
#' @param min_size Integer. Minimum number of genes per gene set.
#' @param max_size Integer. Maximum number of genes per gene set.
#' @param gsea_param  Float. GSEA parameter. Defaults to `1.0`.
#' @param sample_size Integer. Number of samples to iterate through for the
#' multi-level implementation of fgsea.
#' @param eps Float. Boundary for calculating the p-value. Used for the multi-
#' level implementation of fgsea.
#'
#' @returns List with parameters for usage in subsequent function.
#'
#' @export
params_gsea <- function(
  min_size = 5L,
  max_size = 500L,
  gsea_param = 1.0,
  sample_size = 101L,
  eps = 1e-50
) {
  # Checks
  checkmate::qassert(min_size, "I1[3,)")
  checkmate::qassert(max_size, "I1[4,)")
  checkmate::qassert(gsea_param, "N1")
  checkmate::qassert(sample_size, "N1")
  checkmate::qassert(eps, "N1")
  # Returns
  return(
    list(
      min_size = min_size,
      max_size = max_size,
      gsea_param = gsea_param,
      sample_size = sample_size,
      eps = eps
    )
  )
}

## gsva ------------------------------------------------------------------------

#' Wrapper function to generate GSVA parameters
#'
#' @param tau Float. Tau parameter, usual recommendation is to use `1.0` here.
#' Larger values emphasise the tails more.
#' @param min_size Integer. Minimum number of genes per gene set.
#' @param max_size Integer. Maximum number of genes per gene set.
#' @param max_diff Boolean. Scoring mode for GSVA, if `TRUE` = difference; if
#' `FALSE` = larger absolute value.
#' @param abs_rank Boolean. If `TRUE` = pos - neg, `FALSE` = pos + neg.
#'
#' @returns List with parameters for usage in subsequent function.
#'
#' @export
params_gsva <- function(
  tau = 1.0,
  min_size = 5L,
  max_size = 500L,
  max_diff = TRUE,
  abs_rank = FALSE
) {
  # checks
  checkmate::qassert(tau, "N1")
  checkmate::qassert(min_size, "I1")
  checkmate::qassert(max_size, "I1")
  checkmate::qassert(max_diff, "B1")
  checkmate::qassert(abs_rank, "B1")

  # return
  return(
    list(
      tau = tau,
      min_size = min_size,
      max_size = max_size,
      max_diff = max_diff,
      abs_rank = abs_rank
    )
  )
}

## ssgsea ----------------------------------------------------------------------

#' Wrapper function to generate ssGSEA parameters
#'
#' @param alpha Float. The exponent defining the weight of the tail in the
#' random walk performed by ssGSEA.
#' @param min_size Integer. Minimum number of genes per gene set.
#' @param max_size Integer. Maximum number of genes per gene set.
#' @param normalise Boolean. Shall the scores be normalised.
#'
#' @returns List with parameters for usage in subsequent function.
#'
#' @export
params_ssgsea <- function(
  alpha = 0.25,
  min_size = 5L,
  max_size = 500L,
  normalise = TRUE
) {
  # checks
  checkmate::qassert(alpha, "N1(0, 1)")
  checkmate::qassert(min_size, "I1")
  checkmate::qassert(max_size, "I1")
  checkmate::qassert(normalise, "B1")

  # return
  return(
    list(
      alpha = alpha,
      min_size = min_size,
      max_size = max_size,
      normalise = TRUE
    )
  )
}

## coremo ----------------------------------------------------------------------

#' Wrapper function to generate CoReMo parameters
#'
#' @param epsilon Float. Epsilon parameter for the chosen RBF function, see
#' `rbf_func`. The higher, the more aggressively low correlations will be
#' shrunk.
#' @param k_min,k_max Integer. Minimum and maximum number of cuts to use for the
#' hierarchical clustering.
#' @param min_size Optional integer. Minimum size of the clusters. Smaller
#' clusters will be combined together.
#' @param junk_module_threshold Float. Threshold for the minimum correlation
#' to be observed in a module. Defaults to `0.05`.
#' @param rbf_func String. Type of RBF you wish to apply to down-weigh weak
#' correlations. Defaults to `"gaussian"`.
#' @param cor_method String. The type of correlation to use. Defaults to
#' `"spearman"`.
#'
#' @returns List with parameters for usage in subsequent function.
#'
#' @export
params_coremo <- function(
  epsilon = 2,
  k_min = 2L,
  k_max = 150L,
  min_size = NULL,
  junk_module_threshold = 0.05,
  rbf_func = c("gaussian", "inverse_quadratic", "bump"),
  cor_method = c("spearman", "pearson")
) {
  # Standard choices
  cor_method <- match.arg(cor_method)
  rbf_func <- match.arg(rbf_func)
  # Checks
  checkmate::qassert(epsilon, "N1")
  checkmate::qassert(k_min, "I1")
  checkmate::qassert(k_max, "I1")
  checkmate::qassert(min_size, c("I1", "0"))
  checkmate::assertChoice(rbf_func, c("gaussian", "inverse_quadratic", "bump"))
  checkmate::assertChoice(
    cor_method,
    c("spearman", "pearson")
  )
  # Returns
  return(
    list(
      epsilon = epsilon,
      k_min = k_min,
      k_max = k_max,
      min_size = min_size,
      junk_module_threshold = junk_module_threshold,
      rbf_func = rbf_func,
      cor_method = cor_method
    )
  )
}

## DGRDL -----------------------------------------------------------------------

#' Wrapper function to generate DGRDL parameters
#'
#' @param sparsity Integer. Sparsity constraint (max non-zero coefficients per
#' signal)
#' @param dict_size Integer. Dictionary size
#' @param alpha Float. Sample context regularisation weight.
#' @param beta Float. Feature effect regularisation weight.
#' @param max_iter Integer. Maximum number of iterations for the main algorithm.
#' @param k_neighbours Integer. Number of neighbours in the KNN graph.
#' @param admm_iter Integer. ADMM iterations for sparse coding.
#' @param rho Float. ADMM step size.
#'
#' @returns List with parameters for usage in subsequent function.
#'
#' @export
params_dgrdl <- function(
  sparsity = 5L,
  dict_size = 5L,
  alpha = 1.0,
  beta = 1.0,
  max_iter = 20L,
  k_neighbours = 5L,
  admm_iter = 5L,
  rho = 1.0
) {
  # checks
  checkmate::qassert(sparsity, "I1")
  checkmate::qassert(dict_size, "I1")
  checkmate::qassert(alpha, "N1")
  checkmate::qassert(beta, "N1")
  checkmate::qassert(max_iter, "I1")
  checkmate::qassert(k_neighbours, "I1")
  checkmate::qassert(admm_iter, "I1")
  checkmate::qassert(rho, "N1")

  list(
    sparsity = sparsity,
    dict_size = dict_size,
    alpha = alpha,
    beta = beta,
    max_iter = max_iter,
    k_neighbours = k_neighbours,
    admm_iter = admm_iter,
    rho = rho
  )
}

## SNF -------------------------------------------------------------------------

#' Wrapper function to generate SNF parameters
#'
#' @param k Integer. Number of neighbours to consider.
#' @param t Integer. Number of iterations for the SNF algorithm.
#' @param mu Float. Normalisation factor for the Gaussian kernel width.
#' @param alpha Float. Normalisation parameter controlling the fusion strength.
#' @param normalise Boolean. Shall continuous values be Z-scored.
#' @param distance_metric String. One of
#' `c("euclidean", "manhattan", "canberra", "cosine")`. Which distance metric
#' to use for the continuous calculations. In case of pure categorical, Hamming
#' will be used, for mixed data types Gower distance is used.
#'
#' @returns List with parameters for usage in subsequent function.
#'
#' @export
params_snf <- function(
  k = 20L,
  t = 20L,
  mu = 0.5,
  alpha = 1.0,
  normalise = TRUE,
  distance_metric = c("euclidean", "manhattan", "canberra", "cosine")
) {
  distance_metric <- match.arg(distance_metric)

  # checks
  checkmate::qassert(k, "I1")
  checkmate::qassert(t, "I1")
  checkmate::qassert(mu, "N1[0, 1]")
  checkmate::qassert(alpha, "N1")
  checkmate::qassert(normalise, "B1")
  checkmate::assertChoice(
    distance_metric,
    c("euclidean", "manhattan", "canberra", "cosine")
  )

  list(
    k = k,
    t = t,
    mu = mu,
    alpha = alpha,
    distance_metric = distance_metric,
    normalise = normalise
  )
}

## single cell -----------------------------------------------------------------

### general --------------------------------------------------------------------

# TODO - too many kNN parameters - refactor them out

### synthetic data -------------------------------------------------------------

#' Default parameters for generation of synthetic data
#'
#' @param n_cells Integer. Number of cells.
#' @param n_genes Integer. Number of genes.
#' @param marker_genes List. A nested list that indicates which gene indices
#' are markers for which cell.
#' @param n_batches Integer. Number of batches.
#' @param batch_effect_strength String. One of
#' `c("strong", "medium", "weak")`. The strength of the batch effect to add.
#'
#' @return A list with the parameters.
#'
#' @export
params_sc_synthetic_data <- function(
  n_cells = 1000L,
  n_genes = 100L,
  marker_genes = list(
    cell_type_1 = list(
      marker_genes = 0:9L
    ),
    cell_type_2 = list(
      marker_genes = 10:19L
    ),
    cell_type_3 = list(
      marker_genes = 20:29L
    )
  ),
  n_batches = 1L,
  batch_effect_strength = c("strong", "medium", "weak")
) {
  batch_effect_strength <- match.arg(batch_effect_strength)

  # checks
  checkmate::qassert(n_cells, "I1")
  checkmate::qassert(n_genes, "I1")
  checkmate::assertList(marker_genes, types = "list", names = "named")
  checkmate::qassert(n_batches, "I1")
  checkmate::assertChoice(batch_effect_strength, c("strong", "medium", "weak"))

  list(
    n_cells = n_cells,
    n_genes = n_genes,
    marker_genes = marker_genes,
    n_batches = n_batches,
    batch_effect_strength = batch_effect_strength
  )
}

### io -------------------------------------------------------------------------

#' Wrapper function to provide data for mtx-based loading
#'
#' @param path_mtx String. Path to the .mtx file
#' @param path_obs String. Path to the file containing cell/barcode info.
#' @param path_var String. Path to the file containing gene/variable info.
#' @param cells_as_rows Boolean. Do cells represent the rows or columns.
#' @param has_hdr Boolean. Do the plain text files have headers.
#'
#' @returns A list with the mtx loading parameters for usage in subsequent
#' functions.
#'
#' @export
params_sc_mtx_io <- function(
  path_mtx,
  path_obs,
  path_var,
  cells_as_rows,
  has_hdr
) {
  # checks
  checkmate::assertFileExists(path_mtx)
  checkmate::assertFileExists(path_obs)
  checkmate::assertFileExists(path_var)
  checkmate::qassert(cells_as_rows, "B1")
  checkmate::qassert(has_hdr, "B1")

  list(
    path_mtx = path_mtx,
    path_obs = path_obs,
    path_var = path_var,
    cells_as_rows = cells_as_rows,
    has_hdr = has_hdr
  )
}

### qc -------------------------------------------------------------------------

#' Wrapper function to generate QC metric params for single cell
#'
#' @param min_unique_genes Integer. Minimum number of unique genes per cell/spot
#' to be included.
#' @param min_lib_size Integer. Minimum library size per cell/spot to be
#' included.
#' @param min_cells Integer. Minimum number of cells a gene has to be
#' expressed to be included.
#' @param target_size Float. The target size for the normalisation. Defaults
#' to `1e4`.
#'
#' @returns A list with the minimum quality parameters + target size.
#'
#' @export
params_sc_min_quality <- function(
  min_unique_genes = 100L,
  min_lib_size = 250L,
  min_cells = 10L,
  target_size = 1e4
) {
  # checks
  checkmate::qassert(min_unique_genes, "I1")
  checkmate::qassert(min_lib_size, "I1")
  checkmate::qassert(min_cells, "I1")
  checkmate::qassert(target_size, "N1")

  list(
    min_unique_genes = min_unique_genes,
    min_lib_size = min_lib_size,
    min_cells = min_cells,
    target_size = target_size
  )
}

#' Wrapper function for HVG detection parameters.
#'
#' @param method String. One of `c("vst", "meanvarbin", "dispersion")`.
#' @param loess_span Numeric. The span parameter for the loess function that is
#' used to standardise the variance for `method = "vst"`.
#' @param num_bin Integer. Not yet implemented.
#' @param bin_method String. One of `c("equal_width", "equal_freq")`.
#'
#' @returns A list with the HVG parameters
#'
#' @export
params_sc_hvg <- function(
  method = "vst",
  loess_span = 0.3,
  num_bin = 20L,
  bin_method = "equal_width"
) {
  # check
  checkmate::assertChoice(method, c("vst", "meanvarbin", "dispersion"))
  checkmate::qassert(loess_span, "N1[0.1, 1]")
  checkmate::qassert(num_bin, "N1")
  checkmate::assertChoice(bin_method, c("equal_width", "equal_freq"))

  list(
    method = method,
    loess_span = loess_span,
    num_bin = num_bin,
    bin_method = bin_method
  )
}

#' Wrapper function for parameters for neighbour identification in single cell
#'
#' @param k Integer. Number of neighbours to return.
#' @param knn_algorithm String. One of `c("annoy", "hnsw", "nndescent")`.
#' Defaults to `"annoy"`.
#' @param ann_dist String. One of `c("cosine", "euclidean")`. The distance
#' metric to be used for the approximate neighbour search.
#' @param n_trees Integer. Number of trees to use for the `annoy` algorithm.
#' @param search_budget Integer. Search budget per tree for the `annoy`
#' algorithm.
#' @param max_iter Integer. Maximum iterations for the `nndescent` algorithm.
#' @param rho Numeric. Sampling rate for the `nndescent` algorithm.
#' @param delta Numeric. Early termination criterium for the `nndescent`
#' algorithm.
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
#' lower-ranked (closer) shared neighbors result in higher edge weights
#' Both methods produce weights normalised to the range `[0, 1]`.
#'
#' @returns A list with the neighbour parameters.
#'
#' @export
params_sc_neighbours <- function(
  k = 15L,
  knn_algorithm = c("annoy", "hnsw", "nndescent"),
  ann_dist = c("cosine", "euclidean"),
  n_trees = 100L,
  search_budget = 100L,
  max_iter = 25L,
  rho = 1.0,
  delta = 0.001,
  full_snn = FALSE,
  pruning = 1 / 15,
  snn_similarity = c("rank", "jaccard")
) {
  knn_algorithm <- match.arg(knn_algorithm)
  snn_similarity <- match.arg(snn_similarity)
  ann_dist <- match.arg(ann_dist)

  # check
  checkmate::qassert(k, "I1")
  checkmate::assertChoice(knn_algorithm, c("annoy", "hnsw", "nndescent"))
  checkmate::assertChoice(ann_dist, c("cosine", "euclidean"))
  checkmate::qassert(n_trees, "I1")
  checkmate::qassert(search_budget, "N1")
  checkmate::qassert(max_iter, "I1")
  checkmate::qassert(rho, "N1")
  checkmate::qassert(delta, "N1")
  checkmate::qassert(full_snn, "B1")
  checkmate::qassert(pruning, "N1[0, 1]")
  checkmate::assertChoice(snn_similarity, c("rank", "jaccard"))

  list(
    k = k,
    knn_algorithm = knn_algorithm,
    ann_dist = ann_dist,
    n_trees = n_trees,
    search_budget = search_budget,
    max_iter = max_iter,
    rho = rho,
    delta = delta,
    full_snn = full_snn,
    pruning = pruning,
    snn_similarity = snn_similarity
  )
}

### batch corrections ----------------------------------------------------------

#### BBKNN ---------------------------------------------------------------------

#' Wrapper function for the BBKNN parameters
#'
#' @param neighbours_within_batch Integer. Number of neighbours to consider
#' per batch.
#' @param knn_method String. One of `c("annoy", "hnsw")`. Defaults to
#' `"annoy"`.
#' @param ann_dist String. One of `c("cosine", "euclidean")`. The distance
#' metric to be used for the approximate neighbour search.
#' @param set_op_mix_ratio Numeric. Mixing ratio between union (1.0) and
#' intersection (0.0).
#' @param local_connectivity Numeric. UMAP connectivity computation parameter,
#' how many nearest neighbours of each cell are assumed to be fully connected.
#' @param annoy_n_trees Integer. Number of trees to use in the generation of
#' the Annoy index.
#' @param search_budget Integer. Search budget per tree for the `annoy`
#' algorithm.
#' @param trim Optional integer. Trim the neighbours of each cell to these many
#' top connectivities. May help with population independence and improve the
#' tidiness of clustering. If `NULL`, it defaults to
#' `10 * neighbours_within_batch`.
#'
#' @returns A list with the BBKNN parameters.
#'
#' @export
params_sc_bbknn <- function(
  neighbours_within_batch = 5L,
  knn_method = c("annoy", "hnsw"),
  ann_dist = c("cosine", "euclidean"),
  set_op_mix_ratio = 1.0,
  local_connectivity = 1.0,
  annoy_n_trees = 100L,
  search_budget = 100L,
  trim = NULL
) {
  knn_method <- match.arg(knn_method)
  ann_dist <- match.arg(ann_dist)

  # checks
  checkmate::qassert(neighbours_within_batch, "I1")
  checkmate::assertChoice(knn_method, c("hnsw", "annoy"))
  checkmate::assertChoice(ann_dist, c("cosine", "euclidean"))
  checkmate::qassert(set_op_mix_ratio, "N1[0, 1]")
  checkmate::qassert(local_connectivity, "N1")
  checkmate::qassert(annoy_n_trees, "I1")
  checkmate::qassert(search_budget, "I1")
  checkmate::qassert(trim, c("0", "I1"))

  # returns
  list(
    neighbours_within_batch = neighbours_within_batch,
    knn_method = knn_method,
    ann_dist = ann_dist,
    set_op_mix_ratio = set_op_mix_ratio,
    local_connectivity = local_connectivity,
    annoy_n_trees = annoy_n_trees,
    search_budget = search_budget,
    trim = trim
  )
}

#### fastMNN -------------------------------------------------------------------

#' Wrapper function for the fastMNN parameters
#'
#' @param k Integer. Number of mutual nearest neighbours to identify.
#' Defaults to `20L`.
#' @param sigma Numeric. Bandwidth of the Gaussian smoothing kernel
#' (as proportion of space radius). Defaults to `0.1`.
#' @param knn_method String. One of `c("annoy", "hnsw")`. Defaults to
#' `"annoy"`.
#' @param dist_metric String. One of `c("cosine", "euclidean")`. Defaults to
#' `"cosine"`.
#' @param annoy_n_trees Integer. Number of trees for Annoy index. Defaults to
#' `100L`.
#' @param annoy_search_budget Integer. Search budget per tree for Annoy.
#' Defaults to `100L`.
#' @param cos_norm Logical. Apply cosine normalisation before computing
#' distances. Defaults to `TRUE`.
#' @param var_adj Logical. Apply variance adjustment to avoid kissing effects.
#' Defaults to `TRUE`.
#' @param no_pcs Integer. Number of PCs to use for MNN calculations.
#' Defaults to `30L`.
#' @param random_svd Logical. Use randomised SVD. Defaults to `TRUE`.
#'
#' @returns A list with the fastMNN parameters.
#'
#' @export
params_sc_fastmnn <- function(
  k = 20L,
  sigma = 0.1,
  knn_method = c("annoy", "hnsw"),
  dist_metric = c("cosine", "euclidean"),
  annoy_n_trees = 100L,
  annoy_search_budget = 100L,
  cos_norm = TRUE,
  var_adj = TRUE,
  no_pcs = 30L,
  random_svd = TRUE
) {
  knn_method <- match.arg(knn_method)
  dist_metric <- match.arg(dist_metric)

  checkmate::qassert(k, "I1")
  checkmate::qassert(sigma, "N1")
  checkmate::assertChoice(knn_method, c("annoy", "hnsw"))
  checkmate::assertChoice(dist_metric, c("cosine", "euclidean"))
  checkmate::qassert(annoy_n_trees, "I1")
  checkmate::qassert(annoy_search_budget, "I1")
  checkmate::qassert(cos_norm, "B1")
  checkmate::qassert(var_adj, "B1")
  checkmate::qassert(no_pcs, "I1")
  checkmate::qassert(random_svd, "B1")

  list(
    k = k,
    sigma = sigma,
    knn_method = knn_method,
    dist_metric = dist_metric,
    annoy_n_trees = annoy_n_trees,
    annoy_search_budget = annoy_search_budget,
    cos_norm = cos_norm,
    var_adj = var_adj,
    no_pcs = no_pcs,
    random_svd = random_svd
  )
}

### meta cells -----------------------------------------------------------------

#### meta cell (hdWGCNA) -------------------------------------------------------

#' Wrapper function for parameters for meta cell generation
#'
#' @param max_shared Integer. Maximum number of allowed shared neighbours for
#' the meta cell to be considered. Defaults to `15L`.
#' @param target_no_metacells Integer. Target number of meta-cells to generate.
#' Defaults to `1000L`.
#' @param max_iter Integer. Maximum number of iterations for the algorithm.
#' Defaults to `5000L`.
#' @param k Integer. Number of neighbours to return. Defaults to `25L`.
#' @param knn_method String. One of `c("annoy", "hnsw", "nndescent")`. Defaults
#' to `"annoy"`.
#' @param ann_dist String. One of `c("cosine", "euclidean")`. The distance
#' metric to be used for the approximate neighbour search. Defaults to
#' `"cosine"`.
#' @param n_trees Integer. Number of trees to use for the `annoy` algorithm.
#' Defaults to `100L`.
#' @param search_budget Integer. Search budget per tree for the `annoy`
#' algorithm. Defaults to `100L`.
#' @param nn_max_iter Integer. Maximum iterations for NN Descent. Defaults to
#' `15L`.
#' @param rho Numeric. Sampling rate for NN Descent. Defaults to `1.0`.
#' @param delta Numeric. Early termination criterion for NN Descent. Defaults to
#' `0.001`.
#'
#' @returns A list with the metacell parameters.
#'
#' @export
params_sc_metacells <- function(
  max_shared = 15L,
  target_no_metacells = 1000L,
  max_iter = 5000L,
  k = 25L,
  knn_method = c("annoy", "hnsw", "nndescent"),
  ann_dist = c("cosine", "euclidean"),
  n_trees = 100L,
  search_budget = 100L,
  nn_max_iter = 15L,
  rho = 1.0,
  delta = 0.001
) {
  knn_method <- match.arg(knn_method)
  ann_dist <- match.arg(ann_dist)

  checkmate::qassert(max_shared, "I1")
  checkmate::qassert(target_no_metacells, "I1")
  checkmate::qassert(max_iter, "I1")
  checkmate::qassert(k, "I1")
  checkmate::assertChoice(knn_method, c("annoy", "hnsw", "nndescent"))
  checkmate::assertChoice(ann_dist, c("cosine", "euclidean"))
  checkmate::qassert(n_trees, "I1")
  checkmate::qassert(search_budget, "I1")
  checkmate::qassert(nn_max_iter, "I1")
  checkmate::qassert(rho, "N1")
  checkmate::qassert(delta, "N1")

  list(
    max_shared = max_shared,
    target_no_metacells = target_no_metacells,
    max_iter = max_iter,
    k = k,
    knn_method = knn_method,
    ann_dist = ann_dist,
    n_trees = n_trees,
    search_budget = search_budget,
    nn_max_iter = nn_max_iter,
    rho = rho,
    delta = delta
  )
}

#### sea cells -----------------------------------------------------------------

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
#' @param k Integer. Number of neighbours for the kNN algorithm. Defaults to
#' `25L`.
#' @param knn_method String. One of `c("annoy", "hnsw", "nndescent")`. Defaults
#' to `"annoy"`.
#' @param n_trees Integer. Number of trees for Annoy index. Defaults to `100L`.
#' @param search_budget Integer. Search budget during querying. Defaults to
#' `100L`.
#' @param nn_max_iter Integer. Maximum iterations for NN Descent. Defaults to
#' `15L`.
#' @param rho Numeric. Sampling rate for NN Descent. Defaults to `1.0`.
#' @param delta Numeric. Early termination criterion for NN Descent. Defaults to
#' `0.001`.
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
  k = 25L,
  knn_method = c("annoy", "hnsw", "nndescent"),
  n_trees = 100L,
  search_budget = 100L,
  nn_max_iter = 15L,
  rho = 1.0,
  delta = 0.001
) {
  # checks
  knn_method <- match.arg(knn_method)
  checkmate::qassert(n_sea_cells, "I1")
  checkmate::qassert(max_fw_iters, "I1")
  checkmate::qassert(convergence_epsilon, "N1")
  checkmate::qassert(max_iter, "I1")
  checkmate::qassert(min_iter, "I1")
  checkmate::qassert(greedy_threshold, "I1")
  checkmate::qassert(graph_building, "S1")
  checkmate::qassert(k, "I1")
  checkmate::assertChoice(knn_method, c("annoy", "hnsw", "nndescent"))
  checkmate::qassert(n_trees, "I1")
  checkmate::qassert(search_budget, "I1")
  checkmate::qassert(nn_max_iter, "I1")
  checkmate::qassert(rho, "N1")
  checkmate::qassert(delta, "N1")

  list(
    n_sea_cells = n_sea_cells,
    max_fw_iters = max_fw_iters,
    convergence_epsilon = convergence_epsilon,
    max_iter = max_iter,
    min_iter = min_iter,
    greedy_threshold = greedy_threshold,
    graph_building = graph_building,
    k = k,
    knn_method = knn_method,
    n_trees = n_trees,
    search_budget = search_budget,
    nn_max_iter = nn_max_iter,
    rho = rho,
    delta = delta
  )
}
