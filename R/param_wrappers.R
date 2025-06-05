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
#'
#' @returns List with parameters for usage in subsequent function.
#'
#' @export
params_community_detection <- function(
  max_nodes = 300L,
  min_nodes = 10L,
  min_seed_nodes = 2L,
  initial_res = 0.5
) {
  # Checks
  checkmate::qassert(max_nodes, sprintf("I1[%i,)", min_nodes))
  checkmate::qassert(min_nodes, "I1")
  checkmate::qassert(min_seed_nodes, "I1")
  checkmate::qassert(initial_res, "N1")
  # Return
  return(
    list(
      max_nodes = max_nodes,
      min_nodes = min_nodes,
      min_seed_nodes = min_seed_nodes,
      initial_res = initial_res
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
  sample_size = 121L,
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

## coremo ----------------------------------------------------------------------

#' Wrapper function to generate CoReMo parameters
#'
#' @param k_min,k_max Integer. Minimum and maximum number of cuts to use for the
#' hierarchical clustering.
#' @param min_size Optional integer. Minimum size of the clusters. Smaller
#' clusters will be combined together.
#' @param rbf_func String. Type of RBF you wish to apply to down-weigh weak
#' correlations. Defaults to `"gaussian"`.
#' @param cluster_method String. The type of clustering method you want to use
#' for [stats::hclust()]. Defaults to `"ward.D"`.
#' @param cor_method String. The type of correlation to use. Defaults to
#' `"spearman"`.
#'
#' @returns List with parameters for usage in subsequent function.
#'
#' @export
params_coremo <- function(
  k_min = 2L,
  k_max = 150L,
  min_size = NULL,
  rbf_func = c("gaussian", "inverse_quadratic", "bump"),
  cluster_method = c(
    "ward.D2",
    "ward.D",
    "single",
    "complete",
    "average",
    "mcquitty",
    "median",
    "centroid"
  ),
  cor_method = c("spearman", "pearson")
) {
  # Standard choices
  cluster_method <- match.arg(cluster_method)
  cor_method <- match.arg(cor_method)
  rbf_func <- match.arg(rbf_func)
  # Checks
  checkmate::qassert(k_min, "I1")
  checkmate::qassert(k_max, "I1")
  checkmate::qassert(min_size, c("I1", "0"))
  checkmate::assertChoice(rbf_func, c("gaussian", "inverse_quadratic", "bump"))
  checkmate::assertChoice(
    cluster_method,
    c(
      "ward.D",
      "ward.D2",
      "single",
      "complete",
      "average",
      "mcquitty",
      "median",
      "centroid"
    )
  )
  checkmate::assertChoice(
    cor_method,
    c("spearman", "pearson")
  )
  # Returns
  return(
    list(
      k_min = k_min,
      k_max = k_max,
      min_size = min_size,
      rbf_func = rbf_func,
      cluster_method = cluster_method,
      cor_method = cor_method
    )
  )
}
