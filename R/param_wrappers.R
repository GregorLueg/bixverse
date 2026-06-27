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

## NMF (HALS) ------------------------------------------------------------------

#' Wrapper function for NMF (HALS) parameters
#'
#' @param max_iter Integer. Maximum number of HALS iterations.
#' @param tol Numeric. Convergence tolerance on the relative change in
#' reconstruction loss.
#' @param eps Numeric. Numerical floor for non-negativity / division safety.
#' @param check_every Integer. Convergence check interval in iterations.
#' @param nmf_init String. One of `c("nndsvd", "svd", "random")`. `"nndsvd"`
#' and `"svd"` both map to deterministic NNDSVD initialisation; `"random"`
#' uses random non-negative draws. For stabilised (multi-run) NMF this field
#' is ignored and random init is always used.
#'
#' @returns A list with the HALS NMF parameters.
#'
#' @export
params_nmf_hals <- function(
  max_iter = 250L,
  tol = 1e-4,
  eps = 1e-10,
  check_every = 10L,
  nmf_init = "nndsvd"
) {
  checkmate::qassert(max_iter, "I1[1,)")
  checkmate::qassert(tol, "N1(0,)")
  checkmate::qassert(eps, "N1(0,)")
  checkmate::qassert(check_every, "I1[1,)")
  checkmate::assertChoice(nmf_init, c("nndsvd", "svd", "random"))

  list(
    max_iter = max_iter,
    tol = tol,
    eps = eps,
    check_every = check_every,
    nmf_init = nmf_init
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

## cistarget -------------------------------------------------------------------

#' Wrapper function to CisTarget parameters
#'
#' @param auc_threshold Numeric between 0 and 1. Proportion of genes to use
#' for AUC threshold calculation. Default is 0.05 (5% of genes).
#' @param nes_threshold Numeric. Normalised Enrichment Score threshold for
#' significant motifs. Default is 3.0.
#' @param rcc_method Character. Method for recovery curve calculation. Either
#' "approx" (approximate, faster) or "icistarget" (exact, slower).
#' @param high_conf_cats Character vector. Annotation categories considered
#' high confidence. Default includes direct annotations and orthology-based
#' inferences.
#' @param low_conf_cats Character vector. Annotation categories considered
#' lower confidence. Default includes motif similarity-based inferences.
#'
#' @return A validated list of RcisTarget parameters.
#'
#' @export
params_cistarget <- function(
  auc_threshold = 0.05,
  nes_threshold = 3.0,
  rcc_method = c("approx", "icistarget"),
  high_conf_cats = c("directAnnotation", "inferredBy_Orthology"),
  low_conf_cats = c(
    "inferredBy_MotifSimilarity",
    "inferredBy_MotifSimilarity_n_Orthology"
  )
) {
  rcc_method <- match.arg(rcc_method)

  checkmate::qassert(auc_threshold, "N1[0, 1]")
  checkmate::qassert(nes_threshold, "N1")
  checkmate::assertChoice(rcc_method, c("approx", "icistarget"))
  checkmate::qassert(high_conf_cats, "S+")
  checkmate::qassert(low_conf_cats, "S+")

  return(list(
    auc_threshold = auc_threshold,
    nes_threshold = nes_threshold,
    rcc_method = rcc_method,
    high_conf_cats = high_conf_cats,
    low_conf_cats = low_conf_cats
  ))
}

## graph label propagation -----------------------------------------------------

#' Wrapper function to generate label propagation parameters
#'
#' @param alpha Numeric. Controls the spreading strength. Higher values anchor
#' labelled nodes more strongly to their original label. Defaults to `0.9`.
#' @param iter Integer. Maximum number of iterations to run. Defaults to `100L`.
#' @param tolerance Numeric. Convergence threshold. Stops early if the maximum
#' change across all nodes falls below this value. Defaults to `1e-6`.
#' @param symmetrise Logical. Whether to symmetrise the graph. Defaults to
#' `TRUE`.
#' @param symmetry_strategy Character. Strategy to resolve weight conflicts when
#' symmetrising. One of `"average"`, `"min"`, or `"max"`. Only relevant when
#' `symmetrise = TRUE` and edge weights are provided. Defaults to `"average"`.
#' @param max_hops Integer or NULL. If provided, restricts label spreading to
#' nodes within this many hops of any labelled node. Nodes beyond this limit
#' remain all-zeroes. Defaults to `NULL`.
#'
#' @returns A named list of label propagation parameters.
#'
#' @export
params_label_propagation <- function(
  alpha = 0.9,
  iter = 100L,
  tolerance = 1e-6,
  symmetrise = TRUE,
  symmetry_strategy = "average",
  max_hops = NULL
) {
  checkmate::qassert(alpha, "R1[0, 1]")
  checkmate::qassert(iter, "I1[1,]")
  checkmate::qassert(tolerance, "R1")
  checkmate::qassert(symmetrise, "B1")
  checkmate::assert_choice(symmetry_strategy, c("average", "avg", "min", "max"))
  if (!is.null(max_hops)) {
    checkmate::qassert(max_hops, "I1[0,]")
  }

  list(
    alpha = alpha,
    iter = iter,
    tolerance = tolerance,
    symmetrise = symmetrise,
    symmetry_strategy = symmetry_strategy,
    max_hops = max_hops
  )
}

## single cell -----------------------------------------------------------------

### general --------------------------------------------------------------------

#### synthetic data ------------------------------------------------------------

##### rna ----------------------------------------------------------------------

#' Default parameters for generation of synthetic single cell data (RNA)
#'
#' @description
#' For the generation of synthetic single cell data mostly for testing or
#' showcasing purposes. The default configurations generates 1000 cells x 100
#' genes with genes 1:10 being cell markers for cell type 1, genes 11:20 for
#' cell type 2 and genes 21:30 for cell type.
#'
#' @param n_cells Integer. Number of cells.
#' @param n_genes Integer. Number of genes.
#' @param n_batches Integer. Number of batches.
#' @param marker_genes List. A nested list that indicates which gene indices
#' are markers for which cell.
#' @param batch_effect_strength String. One of
#' `c("strong", "medium", "weak")`. The strength of the batch effect to add.
#' @param n_samples Optional integer. Shall sample membership be added to the
#' synthetic data. If you want sample information you need to provide
#' `n_samples` and `sample_bias`.
#' @param sample_bias Optional string. One of
#' `c("even", "slightly_uneven", "very_uneven")`
#'
#' @return A list with the parameters.
#'
#' @export
params_sc_synthetic_data <- function(
  n_cells = 1000L,
  n_genes = 100L,
  n_batches = 1L,
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
  batch_effect_strength = c("strong", "medium", "weak"),
  n_samples = NULL,
  sample_bias = NULL
) {
  batch_effect_strength <- match.arg(batch_effect_strength)

  # checks
  checkmate::qassert(n_cells, "I1")
  checkmate::qassert(n_genes, "I1")
  checkmate::assertList(marker_genes, types = "list", names = "named")
  checkmate::qassert(n_batches, "I1")
  checkmate::assertChoice(batch_effect_strength, c("strong", "medium", "weak"))
  checkmate::qassert(n_samples, c("I1", "0"))
  checkmate::assert(
    checkmate::testNull(sample_bias),
    checkmate::testChoice(
      sample_bias,
      c("even", "slightly_uneven", "very_uneven")
    )
  )

  list(
    n_cells = n_cells,
    n_genes = n_genes,
    marker_genes = marker_genes,
    n_batches = n_batches,
    batch_effect_strength = batch_effect_strength,
    n_samples = n_samples,
    sample_bias = sample_bias
  )
}

##### adt ----------------------------------------------------------------------

#' Default parameters for generation of synthetic single cell data (ADT)
#'
#' @description
#' For the generation of synthetic single cell data mostly for testing or
#' showcasing purposes. In this case, ADT counts to test multi-modal
#' integration. The default configurations generates 1000 cells x 15
#' proteins with probes 1:3 being cell markers for cell type 1, genes 4:6 for
#' cell type 2 and genes 7:9 for cell type. Columns 13:15 represents isotype
#' controls.
#'
#' @param n_cells Integer. Number of cells.
#' @param n_proteins Integer. Number of proteins
#' @param n_batches Integer. Number of batches.
#' @param marker_genes List. A nested list that indicates which gene indices
#' are markers for which cell.
#' @param isotype_controls Integer vector. The columns that defines the isotype
#' controls. (0-indexed!)
#' @param batch_effect_strength String. One of
#' `c("strong", "medium", "weak")`. The strength of the batch effect to add.
#'
#' @return A list with the parameters.
#'
#' @export
#'
#' @keywords internal
params_sc_synthetic_data_adt <- function(
  n_cells = 1000L,
  n_proteins = 15L,
  n_batches = 1L,
  marker_genes = list(
    cell_type_1 = list(
      marker_genes = 0:2L
    ),
    cell_type_2 = list(
      marker_genes = 3:5L
    ),
    cell_type_3 = list(
      marker_genes = 6:8L
    )
  ),
  isotype_controls = 12L:14L,
  batch_effect_strength = c("strong", "medium", "weak")
) {
  batch_effect_strength <- match.arg(batch_effect_strength)

  # checks
  checkmate::qassert(n_cells, "I1")
  checkmate::qassert(n_proteins, "I1")
  checkmate::assertList(marker_genes, types = "list", names = "named")
  checkmate::qassert(n_batches, "I1")
  checkmate::qassert(isotype_controls, "I+")
  checkmate::assertChoice(batch_effect_strength, c("strong", "medium", "weak"))

  list(
    n_cells = n_cells,
    n_proteins = n_proteins,
    marker_genes = marker_genes,
    n_batches = n_batches,
    isotype_controls = isotype_controls,
    batch_effect_strength = batch_effect_strength
  )
}

#### io ------------------------------------------------------------------------

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
    path_mtx = path.expand(path_mtx),
    path_obs = path.expand(path_obs),
    path_var = path.expand(path_var),
    cells_as_rows = cells_as_rows,
    has_hdr = has_hdr
  )
}

#### qc ------------------------------------------------------------------------

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

#### hvg -----------------------------------------------------------------------

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

#### pca -----------------------------------------------------------------------

#' Wrapper for PCA specifically designed for single cells
#'
#' @param mean_center Boolean. Shall the data be mean centered
#' @param normalise_variance Boolean. Shall the data have normalised variance
#' @param randomised Boolean. Shall fast, approximate randomised SVD be used.
#' @param clr Boolean. Shall the CLR-type `PFlogPF` be applied, see Booeshaghi,
#' et al.
#' @param size_factor Numeric. The used size factor during I/O. It needs to be
#' the same as during I/O to have correct results when using the `PFlogPF`
#' transformation.
#'
#' @returns A list with the parameters
#'
#' @export
params_sc_pca <- function(
  mean_center = TRUE,
  normalise_variance = TRUE,
  randomised = TRUE,
  clr = FALSE,
  size_factor = 1e4
) {
  # checks
  checkmate::qassert(mean_center, "B1")
  checkmate::qassert(normalise_variance, "B1")
  checkmate::qassert(randomised, "B1")
  checkmate::qassert(clr, "B1")
  checkmate::qassert(size_factor, "N1")

  list(
    mean_center = mean_center,
    normalise_variance = normalise_variance,
    randomised = randomised,
    clr = clr,
    size_factor = size_factor
  )
}

#### knn -----------------------------------------------------------------------

#' Parameters for single cell kNN searches
#'
#' @param k Integer. Number of neighbours. Defaults to `15L`.
#' @param knn_method String. Which method to use for the approximate nearest
#' neighbour search. Defaults to `"kmknn"`. One of
#' `c("kmknn", "hnsw", "annoy", "nndescent", "ivf", "exhaustive")`.
#' @param ann_dist String. Distance metric to use. Defaults to `"euclidean"`.
#' One of `c("cosine", "euclidean")`.
#' @param n_trees Integer. Annoy param: number of trees. Defaults to `50L`.
#' @param search_budget Integer or `NULL`. Annoy param: optional search budget
#' per tree. If `NULL`, defaults to `n_trees * k * 20L` internally.
#' @param delta Numeric. NNDescent param: early termination criterion.
#' Defaults to `0.001`.
#' @param diversify_prob Numeric. NNDescent param: diversification probability
#' applied at the end of index construction. Defaults to `0.0`.
#' @param ef_budget Integer or `NULL`. NNDescent param: optional query budget.
#' Higher values improve recall at the cost of speed.
#' @param m Integer. HNSW param: number of connections between layers.
#' Defaults to `16L`.
#' @param ef_construction Integer. HNSW param: size of the dynamic candidate
#' list during construction. Defaults to `200L`.
#' @param ef_search Integer. HNSW param: size of the candidate list at query
#' time. Higher values improve recall at the cost of speed. Defaults to `100L`.
#' @param n_list Integer or `NULL`. IVF param: number of clusters to generate.
#' If `NULL`, defaults to `sqrt(n)` internally.
#' @param n_probe Integer or `NULL`. IVF param: number of clusters to query.
#' If `NULL`, defaults to `sqrt(n_list)` internally.
#'
#' @return A list with the kNN parameters.
#'
#' @export
params_sc_knn <- function(
  k = 15L,
  knn_method = "kmknn",
  ann_dist = "euclidean",
  n_trees = 50L,
  search_budget = NULL,
  delta = 0.001,
  diversify_prob = 0.0,
  ef_budget = NULL,
  m = 16L,
  ef_construction = 200L,
  ef_search = 100L,
  n_list = NULL,
  n_probe = NULL
) {
  checkmate::qassert(k, "I1[1,)")
  checkmate::assertChoice(
    knn_method,
    c("kmknn", "hnsw", "annoy", "nndescent", "ivf", "exhaustive")
  )
  checkmate::assertChoice(ann_dist, c("cosine", "euclidean"))
  checkmate::qassert(n_trees, "I1[1,)")
  checkmate::assert(
    checkmate::checkNull(search_budget),
    checkmate::checkInt(search_budget, lower = 1L)
  )
  checkmate::qassert(delta, "N1(0,)")
  checkmate::qassert(diversify_prob, "N1[0,1]")
  checkmate::assert(
    checkmate::checkNull(ef_budget),
    checkmate::checkInt(ef_budget, lower = 1L)
  )
  checkmate::qassert(m, "I1[1,)")
  checkmate::qassert(ef_construction, "I1[1,)")
  checkmate::qassert(ef_search, "I1[1,)")
  checkmate::assert(
    checkmate::checkNull(n_list),
    checkmate::checkInt(n_list, lower = 1L)
  )
  checkmate::assert(
    checkmate::checkNull(n_probe),
    checkmate::checkInt(n_probe, lower = 1L)
  )

  res <- list(
    k = k,
    knn_method = knn_method,
    ann_dist = ann_dist,
    n_trees = n_trees,
    search_budget = search_budget,
    delta = delta,
    diversify_prob = diversify_prob,
    ef_budget = ef_budget,
    m = m,
    ef_construction = ef_construction,
    ef_search = ef_search,
    n_list = n_list,
    n_probe = n_probe
  )

  res
}

## single cell (multi modal) ---------------------------------------------------

### dsb adt normalisation ------------------------------------------------------

#' Default parameters for DSB ADT normalisation
#'
#' @param denoise_counts Boolean. Run Step II (cell-to-cell technical noise
#' removal).
#' @param use_isotype_controls Boolean. Include isotype controls in the noise
#' matrix in Step II. Requires `isotype_indices` to be passed at call time.
#' @param pseudocount Numeric. Pseudocount added before the log transform.
#' The DSB paper recommends `10` with empty droplets and `1` without.
#' @param quantile_low Optional numeric in `[0, 1)`. Lower quantile for
#' per-protein output clipping. If `NULL` (and `quantile_high` is also `NULL`),
#' no clipping is applied.
#' @param quantile_high Optional numeric in `(0, 1]`. Upper quantile for
#' per-protein output clipping. If `NULL` (and `quantile_low` is also `NULL`),
#' no clipping is applied.
#'
#' @return A list with the parameters.
#'
#' @export
params_sc_dsb <- function(
  denoise_counts = TRUE,
  use_isotype_controls = TRUE,
  pseudocount = 10,
  quantile_low = NULL,
  quantile_high = NULL
) {
  # checks
  checkmate::qassert(denoise_counts, "B1")
  checkmate::qassert(use_isotype_controls, "B1")
  checkmate::qassert(pseudocount, "N1(0,)")
  checkmate::qassert(quantile_low, c("N1[0,1)", "0"))
  checkmate::qassert(quantile_high, c("N1(0,1]", "0"))

  # both-or-neither for clipping
  if (xor(is.null(quantile_low), is.null(quantile_high))) {
    stop("quantile_low and quantile_high must both be provided or both NULL.")
  }
  if (
    !is.null(quantile_low) &&
      !is.null(quantile_high) &&
      quantile_low >= quantile_high
  ) {
    stop("quantile_low must be strictly less than quantile_high.")
  }

  list(
    denoise_counts = denoise_counts,
    use_isotype_controls = use_isotype_controls,
    pseudocount = pseudocount,
    quantile_low = quantile_low,
    quantile_high = quantile_high
  )
}

### symphony -------------------------------------------------------------------

#' Default parameters for Symphony query mapping
#'
#' @param sigma Numeric. Soft-clustering fuzziness for query -> reference
#' centroid assignment. Symphony R default is 0.1.
#' @param lambda Numeric. Ridge penalty on batch coefficients. Symphony R
#' hardcodes 1.0.
#'
#' @return A list with the parameters.
#'
#' @export
params_symphony_map <- function(sigma = 0.1, lambda = 1.0) {
  checkmate::qassert(sigma, "N1[0,)")
  checkmate::qassert(lambda, "N1[0,)")
  res <- list(sigma = sigma, lambda = lambda)
  class(res) <- c("params_symphony_map", "list")
  res
}

### nichenet -------------------------------------------------------------------

#' Parameters for ligand to target influence computation
#'
#' @param lr_sig_hub Numeric in `[0, 1]`. Hub correction strength for the
#' ligand-receptor / signalling layer. 0 disables correction.
#' @param gr_hub Numeric in `[0, 1]`. Hub correction strength for the gene
#' regulatory layer. 0 disables correction.
#' @param ltf_cutoff Numeric in `[0, 1]`. Quantile cutoff applied to the
#' intermediate ligand-to-TF matrix.
#' @param damping_factor Numeric in `[0, 1]`. PageRank-style damping factor.
#' @param tol Numeric > 0. Convergence tolerance for the propagation step.
#' @param max_iter Integer >= 1. Maximum iterations for the propagation step.
#' @param topology_correction Boolean. Apply topology correction.
#' @param secondary_targets Boolean. Run a second round through targets.
#'
#' @returns A named list of parameters.
#'
#' @export
params_ligand_target <- function(
  lr_sig_hub = 0,
  gr_hub = 0,
  ltf_cutoff = 0.99,
  damping_factor = 0.5,
  tol = 1e-6,
  max_iter = 1000L,
  topology_correction = FALSE,
  secondary_targets = FALSE
) {
  checkmate::qassert(lr_sig_hub, "N1[0,1]")
  checkmate::qassert(gr_hub, "N1[0,1]")
  checkmate::qassert(ltf_cutoff, "N1[0,1]")
  checkmate::qassert(damping_factor, "N1[0,1]")
  checkmate::qassert(tol, "N1(0,)")
  checkmate::qassert(max_iter, "X1[1,)")
  checkmate::qassert(topology_correction, "B1")
  checkmate::qassert(secondary_targets, "B1")

  list(
    lr_sig_hub = lr_sig_hub,
    gr_hub = gr_hub,
    ltf_cutoff = ltf_cutoff,
    damping_factor = damping_factor,
    tol = tol,
    max_iter = as.integer(max_iter),
    topology_correction = topology_correction,
    secondary_targets = secondary_targets
  )
}
