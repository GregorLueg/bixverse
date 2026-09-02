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

## NMF (consensus) -------------------------------------------------------------

#' Wrapper function for consensus NMF parameters
#'
#' @description
#' Controls the clustering step of consensus NMF: the components of every
#' restart are pooled, outliers are dropped by local density, and the survivors
#' are k-means clustered into `k` groups whose median becomes the consensus
#' factor.
#'
#' @param consensus_target String. One of `c("h", "w")`. `"h"` clusters the
#' gene programmes (the spectra), which is what cNMF does and what you almost
#' always want. `"w"` clusters in sample space instead, which on single cell
#' data means cell space: it pools a dense `(k * n_runs) x n_cells` matrix and
#' runs an exhaustive cosine search over it, so it gets expensive fast.
#' @param n_neighbours Integer. Neighbours used for the local density estimate.
#' `0L` picks `ceiling(0.3 * n_runs)` for you.
#' @param density_threshold Numeric. Components whose mean cosine distance to
#' their neighbours exceeds this are dropped as unstable. Cosine distance
#' cannot exceed 2, so any value `>= 2` disables the filter entirely.
#' @param kmeans_iters Integer. Maximum k-means iterations.
#' @param kmeans_n_init Integer. Number of k-means restarts.
#'
#' @returns A list with the consensus NMF parameters.
#'
#' @references Kotliar et al., eLife, 2019
#'
#' @export
params_nmf_consensus <- function(
  consensus_target = c("h", "w"),
  n_neighbours = 0L,
  density_threshold = 0.5,
  kmeans_iters = 100L,
  kmeans_n_init = 3L
) {
  consensus_target <- match.arg(consensus_target)

  checkmate::assertChoice(consensus_target, c("h", "w"))
  checkmate::qassert(n_neighbours, "I1[0,)")
  checkmate::qassert(density_threshold, "N1[0,2]")
  checkmate::qassert(kmeans_iters, "I1[1,)")
  checkmate::qassert(kmeans_n_init, "I1[1,)")

  list(
    consensus_target = consensus_target,
    n_neighbours = n_neighbours,
    density_threshold = density_threshold,
    kmeans_iters = kmeans_iters,
    kmeans_n_init = kmeans_n_init
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
#' @description
#' `auc_threshold` and `max_rank` are two different cutoffs and are easy to
#' confuse. The first truncates the recovery curve used to score motif
#' enrichment, the second sets how deep into the ranking the background curve
#' (mean + 2 SD across all motifs) is built, and therefore where the leading
#' edge cuts. RcisTarget uses `maxRank = 5000` and `nMean = 100`; the defaults
#' here follow it. `auc_threshold` follows pySCENIC at 5%, RcisTarget itself
#' uses 3%.
#'
#' @param auc_threshold Numeric between 0 and 1. Proportion of genes to use
#' for AUC threshold calculation. Default is 0.05 (5% of genes).
#' @param nes_threshold Numeric. Normalised Enrichment Score threshold for
#' significant motifs. Default is 3.0.
#' @param max_rank Integer. Depth of the recovery curves used to derive the
#' background and the leading edge. Clamped to the number of genes in the
#' ranking database. Default is 5000, the RcisTarget value.
#' @param n_mean Integer. Window for the rolling mean smoothing the background
#' recovery curve. Only read when `rcc_method = "approx"`. Default is 100, the
#' RcisTarget value.
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
  max_rank = 5000L,
  n_mean = 100L,
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
  checkmate::qassert(max_rank, "I1[1,)")
  checkmate::qassert(n_mean, "I1[1,)")
  checkmate::assertChoice(rcc_method, c("approx", "icistarget"))
  checkmate::qassert(high_conf_cats, "S+")
  checkmate::qassert(low_conf_cats, "S+")

  return(list(
    auc_threshold = auc_threshold,
    nes_threshold = nes_threshold,
    max_rank = max_rank,
    n_mean = n_mean,
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

## module membership (matrix factorisations) ------------------------------------

#' Wrapper function to generate module membership parameters
#'
#' @description
#' Controls how a `gene x k` loading matrix from ICA, NMF or DGRDL is turned into
#' module membership. Genes are kept where they sit in the tail of a component's
#' loading distribution, which means membership is not exclusive: a gene loading
#' strongly on three components belongs to three modules. Genes in no tail belong
#' to nothing, which is the background category an argmax assignment cannot give
#' you.
#'
#' @param method String. `"zscore"` standardises each component and keeps
#' `abs(z) > cutoff`. `"fdr"` converts to two-sided p-values against a Normal
#' null fitted the same way, Benjamini-Hochberg adjusts, and keeps
#' `padj < fdr`. Defaults to `"zscore"`.
#' @param cutoff Float. Absolute z threshold for `method = "zscore"`. Defaults
#' to `3.0`.
#' @param fdr Float. Adjusted p-value threshold for `method = "fdr"`. Defaults
#' to `0.05`.
#' @param tails String. `"auto"` uses an upper-tail-only test when every loading
#' is non-negative (the NMF case) and a two-sided one otherwise. `"upper"` and
#' `"both"` force the choice. Defaults to `"auto"`.
#' @param scaling String. `"robust"` centres and scales each component by its
#' median and MAD. `"standard"` uses the mean and standard deviation instead,
#' which is stricter and less forgiving of skewed loadings (e.g. NMF).
#' Defaults to `"robust"`.
#'
#' @returns A list with the parameters for usage in subsequent functions.
#'
#' @references Biton, et al., Cell Rep, 2014
#'
#' @export
params_module_membership <- function(
  method = c("zscore", "fdr"),
  cutoff = 3.0,
  fdr = 0.05,
  tails = c("auto", "upper", "both"),
  scaling = c("robust", "standard")
) {
  # Standard choices
  method <- match.arg(method)
  tails <- match.arg(tails)
  scaling <- match.arg(scaling)
  # Checks
  checkmate::assertChoice(method, c("zscore", "fdr"))
  checkmate::assertChoice(tails, c("auto", "upper", "both"))
  checkmate::assertChoice(scaling, c("robust", "standard"))
  checkmate::qassert(cutoff, "N1(0,)")
  checkmate::qassert(fdr, "N1(0,1]")
  # Return
  return(
    list(
      method = method,
      cutoff = cutoff,
      fdr = fdr,
      tails = tails,
      scaling = scaling
    )
  )
}

## synthetic data (bulk) -------------------------------------------------------

### generator ------------------------------------------------------------------

#' Wrapper function to generate synthetic bulk RNAseq parameters
#'
#' @description
#' Parameters for [bixverse::synthetic_bulk_cor_matrix()]. Counts come from a
#' negative binomial with a mean-dispersion trend; co-expression modules are
#' planted by putting each module's genes on a shared latent factor. The
#' `generator` picks how loadings and factors are drawn, which is what makes a
#' given dataset a fair or unfair benchmark for a given method:
#'
#' \itemize{
#'  \item `"hub_modular"` - LogNormal loadings on a Normal factor. Some genes
#'  end up far more connected than others, so this is the WGCNA-style default.
#'  \item `"modular"` - Beta(5, 2) loadings on a Normal factor. Homogeneous
#'  within-module correlation and no hubs.
#'  \item `"non_negative_factor"` - LogNormal loadings on a Gamma factor. The
#'  activity matrix is non-negative by construction, so NMF has a ground truth
#'  it can actually reach.
#'  \item `"non_gaussian_factor"` - LogNormal loadings on a Laplace factor.
#'  Non-Gaussian sources satisfy ICA identifiability.
#' }
#'
#' @param num_samples Integer. Number of samples (columns) to simulate.
#' @param num_genes Integer. Number of genes (rows) to simulate.
#' @param module_sizes Integer vector. Sizes of the co-expression modules. The
#' sum must be smaller or equal to `num_genes`. Genes are assigned in
#' contiguous blocks from the first gene onwards, any remainder is background.
#' Use `integer(0)` for no modules. Must be an integer vector, see the note
#' below.
#' @param generator String. Which topology and distribution family to plant.
#' One of `c("hub_modular", "modular", "non_negative_factor",
#' "non_gaussian_factor")`, see the description. Defaults to `"hub_modular"`.
#' @param seed Integer. Seed for reproducibility purposes.
#' @param mean_exp_gamma_shape,mean_exp_gamma_scale Float. Shape and scale of
#' the Gamma the per-gene mean expression is drawn from.
#' @param disp_intercept,disp_slope Float. Intercept and slope of the negative
#' binomial dispersion trend `disp = 1 / (a + b * mean)`. This is what gives you
#' heteroskedasticity: lowly expressed genes show higher variance.
#' @param noise_std Float. Per-gene per-sample noise standard deviation on the
#' latent log-signal. Smaller values track the module factor more tightly and
#' give stronger within-module correlation. Defaults to `0.1` rather than the
#' crate's `0.3`, see the note below.
#' @param factor_std Float. Standard deviation of the Normal factor. Only used
#' by `"hub_modular"` and `"modular"`; the other two generators draw their
#' factor from `factor_shape`/`factor_scale` instead. Defaults to `0.5` rather
#' than the crate's `0.3`, see the note below.
#' @param factor_shape,factor_scale Float. Shape and scale of the Gamma factor
#' for `"non_negative_factor"`. `factor_scale` doubles as the Laplace scale for
#' `"non_gaussian_factor"`.
#' @param loading_mu,loading_sigma Float. Location and scale of the LogNormal
#' the loadings are drawn from. Unused by `"modular"`, which draws Beta(5, 2).
#' @param hub_percentile Float. Top fraction of module genes flagged as hubs by
#' loading rank. Must be in `(0, 1]`.
#'
#' @returns A list with the parameters for usage in subsequent functions.
#'
#' @details
#' `noise_std` and `factor_std` default to `0.1` and `0.5`, not to the
#' `bixverse-rs` values of `0.3` and `0.3`. At the crate defaults the
#' `"modular"` generator plants modules too weakly to detect at 1000 genes by
#' 100 samples: the within-module minus cross-module mean absolute Spearman gap
#' comes out around `0.06`, against `0.17` to `0.23` for the other three. The
#' values here put all four generators in the `0.30` to `0.39` band, so a
#' comparison across generators reflects the method rather than the signal
#' strength it happened to be handed. Pass the crate values explicitly if you
#' want a harder problem.
#'
#' @references Zhang & Horvath, Stat Appl Genet Mol Biol, 2005
#'
#' @export
params_synthetic_bulk_rnaseq <- function(
  num_samples = 100L,
  num_genes = 1000L,
  module_sizes = c(100L, 100L, 100L),
  generator = c(
    "hub_modular",
    "modular",
    "non_negative_factor",
    "non_gaussian_factor"
  ),
  seed = 123L,
  mean_exp_gamma_shape = 5.0,
  mean_exp_gamma_scale = 10.0,
  disp_intercept = 0.2,
  disp_slope = 0.3,
  noise_std = 0.1,
  factor_std = 0.5,
  factor_shape = 2.0,
  factor_scale = 0.3,
  loading_mu = 0.0,
  loading_sigma = 0.7,
  hub_percentile = 0.1
) {
  # Standard choices
  generator <- match.arg(generator)
  # Checks
  checkmate::qassert(num_samples, "I1[1,)")
  checkmate::qassert(num_genes, "I1[1,)")
  # `module_sizes` reaches Rust via `as_integer_slice()`. A double vector reads
  # as absent there and silently falls back to the crate default, so the
  # integer type is load-bearing rather than cosmetic. Same story for
  # `generator`, where an unknown string quietly resolves to "hub_modular".
  checkmate::qassert(module_sizes, "I*")
  checkmate::assertChoice(
    generator,
    c("hub_modular", "modular", "non_negative_factor", "non_gaussian_factor")
  )
  checkmate::qassert(seed, "I1[0,)")
  checkmate::qassert(mean_exp_gamma_shape, "N1(0,)")
  checkmate::qassert(mean_exp_gamma_scale, "N1(0,)")
  checkmate::qassert(disp_intercept, "N1(0,)")
  checkmate::qassert(disp_slope, "N1(0,)")
  checkmate::qassert(noise_std, "N1[0,)")
  checkmate::qassert(factor_std, "N1(0,)")
  checkmate::qassert(factor_shape, "N1(0,)")
  checkmate::qassert(factor_scale, "N1(0,)")
  checkmate::qassert(loading_mu, "N1")
  checkmate::qassert(loading_sigma, "N1(0,)")
  checkmate::qassert(hub_percentile, "N1(0,1]")
  checkmate::assertTRUE(sum(module_sizes) <= num_genes)
  # Return
  return(
    list(
      num_samples = num_samples,
      num_genes = num_genes,
      module_sizes = module_sizes,
      generator = generator,
      seed = seed,
      mean_exp_gamma_shape = mean_exp_gamma_shape,
      mean_exp_gamma_scale = mean_exp_gamma_scale,
      disp_intercept = disp_intercept,
      disp_slope = disp_slope,
      noise_std = noise_std,
      factor_std = factor_std,
      factor_shape = factor_shape,
      factor_scale = factor_scale,
      loading_mu = loading_mu,
      loading_sigma = loading_sigma,
      hub_percentile = hub_percentile
    )
  )
}

### sparsity -------------------------------------------------------------------

#' Wrapper function to generate bulk sparsification parameters
#'
#' @description
#' Parameters for [bixverse::simulate_dropouts()]. Dropout falls out of the
#' library size rather than an explicit per-gene dropout curve: a size factor
#' `s_j ~ LogNormal(0, capture_efficiency_sigma)` is drawn per sample, giving a
#' target library size of `target_library_size * s_j`, and each gene is
#' binomially thinned towards that target.
#'
#' @param strategy String. Which dropout strategy to apply. Currently only
#' `"seq_depth"`.
#' @param target_library_size Float. Reference library size per sample.
#' @param capture_efficiency_sigma Float. Standard deviation of the LogNormal
#' size-factor distribution. Larger values spread the library sizes further
#' apart.
#' @param seed Integer. Seed for reproducibility purposes.
#'
#' @returns A list with the parameters for usage in subsequent functions.
#'
#' @references Zappia, et al., Genome Biol, 2017
#'
#' @export
params_bulk_sparsity <- function(
  strategy = c("seq_depth"),
  target_library_size = 20000,
  capture_efficiency_sigma = 0.5,
  seed = 123L
) {
  # Standard choices
  strategy <- match.arg(strategy)
  # Checks
  checkmate::assertChoice(strategy, c("seq_depth"))
  checkmate::qassert(target_library_size, "N1(0,)")
  checkmate::qassert(capture_efficiency_sigma, "N1(0,)")
  checkmate::qassert(seed, "I1[0,)")
  # Return
  return(
    list(
      strategy = strategy,
      target_library_size = target_library_size,
      capture_efficiency_sigma = capture_efficiency_sigma,
      seed = seed
    )
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

##### dialogue -----------------------------------------------------------------

#' Default parameters for generation of synthetic DIALOGUE data
#'
#' @description
#' Shapes the fixture with a planted multicellular programme that
#' [bixverse::generate_dialogue_test_data()] builds. The defaults give 14
#' samples x 25 cells x 3 cell types = 1050 cells over 400 genes.
#'
#' @details
#' `n_genes` defaults higher here than on the Rust side, which plants the same
#' structure into 90. R runs the counts through the normal ingestion path, and
#' on a small panel the planted genes are a large enough share of the library
#' that log-normalisation divides the signal back out: at 90 genes the library
#' size tracks the programme at an r of 0.75 and background genes pick up a
#' spurious correlation of their own. At 400 the planted block is a few percent
#' of the library and the contrast survives. The Rust tests feed the planted
#' layer straight in, so they never meet this.
#'
#' @param n_samples Integer. Samples the experiment spans. DIALOGUE needs at
#' least 5.
#' @param cells_per_sample Integer. Cells per sample per cell type.
#' @param n_cell_types Integer. Number of cell types. Must be at least 2.
#' @param n_features Integer. Feature columns per cell type. Must be at least 2.
#' @param n_sample_features Integer. Feature columns carrying a per-sample
#' component. The first of those is the shared programme, the rest are
#' cell-type-specific nuisance; anything past this count is pure noise and
#' exists so the ANOVA filter has something to reject.
#' @param n_genes Integer. Number of genes.
#' @param n_planted Integer. Planted genes per cell type. The blocks are
#' contiguous, so `n_planted * n_cell_types` has to fit into `n_genes`.
#'
#' @return A list with the parameters.
#'
#' @references Jerby-Arnon & Regev, Nature Biotechnology, 2022
#'
#' @export
params_sc_synthetic_dialogue <- function(
  n_samples = 14L,
  cells_per_sample = 25L,
  n_cell_types = 3L,
  n_features = 8L,
  n_sample_features = 5L,
  n_genes = 400L,
  n_planted = 8L
) {
  # checks
  checkmate::qassert(n_samples, "I1[1,)")
  checkmate::qassert(cells_per_sample, "I1[1,)")
  checkmate::qassert(n_cell_types, "I1[2,)")
  checkmate::qassert(n_features, "I1[2,)")
  checkmate::qassert(n_sample_features, "I1[1,)")
  checkmate::qassert(n_genes, "I1[1,)")
  checkmate::qassert(n_planted, "I1[0,)")
  checkmate::assertTRUE(n_sample_features <= n_features)
  checkmate::assertTRUE(n_planted * n_cell_types <= n_genes)

  list(
    n_samples = n_samples,
    cells_per_sample = cells_per_sample,
    n_cell_types = n_cell_types,
    n_features = n_features,
    n_sample_features = n_sample_features,
    n_genes = n_genes,
    n_planted = n_planted
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
#' @param extract_knn Boolean. NNDescent param: hand back the graph the descent
#' already built instead of beam searching it. Skips the query pass entirely,
#' so it is much faster, at the cost of some recall. Rows the descent never
#' filled come back padded with duplicate edges. Ignored by every other method.
#' Defaults to `FALSE`.
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
  extract_knn = FALSE,
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
  checkmate::qassert(extract_knn, "B1")
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
    extract_knn = extract_knn,
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

### sctype ---------------------------------------------------------------------

#' Parameters for the per-cell ScType assignment
#'
#' @description
#' Controls the per-cell path of [assign_sc_type()]: how the raw ScType scores
#' are rescaled, how hard the scores get smoothed over the sNN graph, and where
#' the cut-offs for an Unknown call and for a mixed cluster sit.
#'
#' @param alpha Numeric in `[0, 1]`. Self-retention during smoothing. Each
#' iteration computes `alpha * original + (1 - alpha) * neighbour_average`.
#' @param iterations Integer >= 0. Number of smoothing iterations. `0` disables
#' smoothing.
#' @param tolerance Numeric > 0. Convergence tolerance for the smoothing.
#' @param calibration String. One of `c("none", "column_z")`. `"column_z"`
#' standardises each cell type's score column across cells, which removes the
#' bias towards cell types whose marker sets happen to produce larger scores.
#' @param score_floor Numeric >= 0. Minimum score for a cell to get a call
#' instead of `NA`.
#' @param purity_threshold Numeric in `[0, 1]`. Cluster purity above which the
#' hybrid assignment keeps the cluster-level call.
#'
#' @returns A list with the parameters for usage in subsequent functions.
#'
#' @references Zhou et al., NIPS, 2004.
#'
#' @export
params_sctype_cells <- function(
  alpha = 0.5,
  iterations = 2L,
  tolerance = 1e-4,
  calibration = c("none", "column_z"),
  score_floor = 0.25,
  purity_threshold = 0.9
) {
  calibration <- match.arg(calibration)
  # checks
  checkmate::qassert(alpha, "N1[0,1]")
  checkmate::qassert(iterations, "I1[0,)")
  checkmate::qassert(tolerance, "N1(0,)")
  checkmate::qassert(score_floor, "N1[0,)")
  checkmate::qassert(purity_threshold, "N1[0,1]")
  checkmate::assertChoice(calibration, c("none", "column_z"))
  # return
  res <- list(
    alpha = alpha,
    iterations = as.integer(iterations),
    tolerance = tolerance,
    calibration = calibration,
    score_floor = score_floor,
    purity_threshold = purity_threshold
  )
  class(res) <- c("params_sctype_cells", "list")
  res
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
