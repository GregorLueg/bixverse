# checkmate extensions for various parameter lists provided to the many methods
# in this package

# globals ----------------------------------------------------------------------

# Global to reduce repetition
KNN_PARAM_NAMES <- c(
  "k",
  "knn_method",
  "ann_dist",
  "n_trees",
  "search_budget",
  "delta",
  "diversify_prob",
  "ef_budget",
  "m",
  "ef_construction",
  "ef_search",
  "n_list",
  "n_probe"
)

# internal helpers -------------------------------------------------------------

#' Check that a value is a list with the required names
#'
#' @description Boilerplate guard used at the top of most parameter checkers in
#' this file: verifies `x` is a list and that all `required_names` are present
#' in `names(x)`.
#'
#' @param x The object to check.
#' @param required_names Character vector of names that must be present in
#' `names(x)`.
#'
#' @return `TRUE` if the check was successful, otherwise a checkmate-style error
#' string.
#'
#' @keywords internal
check_list_shape <- function(x, required_names) {
  res <- checkmate::checkList(x)
  if (!isTRUE(res)) {
    return(res)
  }
  checkmate::checkNames(names(x), must.include = required_names)
}

#' Apply qtest rules by name
#'
#' @description Validates the elements of a named list `x` against per-field
#' [checkmate::qtest()] patterns. Fields whose names are not in `rules` are
#' skipped. On failure, returns an error string naming the first offending
#' element and appending an optional `hint`.
#'
#' @param x Named list of parameters to validate.
#' @param rules Named list mapping field name to a qtest pattern (or vector of
#' patterns passed to `qtest`).
#' @param label Short human-readable label used in the error message
#' (e.g. `"GSEA params"`).
#' @param hint Optional string appended to the error message to describe the
#' expected types/ranges. Defaults to `NULL` (no hint).
#'
#' @return `TRUE` if all checked fields pass, otherwise a string of the form
#' `` "The element `<field>` in <label> is invalid. <hint>" ``.
#'
#' @keywords internal
apply_qtest_rules <- function(x, rules, label, hint = NULL) {
  res <- purrr::imap_lgl(x, \(val, name) {
    if (name %in% names(rules)) {
      checkmate::qtest(val, rules[[name]])
    } else {
      TRUE
    }
  })
  if (all(res)) {
    return(TRUE)
  }
  broken <- names(res)[!res][1]
  msg <- sprintf("The element `%s` in %s is invalid.", broken, label)
  if (!is.null(hint)) {
    msg <- paste(msg, hint)
  }
  msg
}

#' Apply testChoice rules by name
#'
#' @description Validates the elements of a named list `x` against per-field
#' [checkmate::testChoice()] choice sets. Fields whose names are not in `rules`
#' are skipped. On failure, returns an error string naming the first offending
#' element and appending an optional `hint`.
#'
#' @param x Named list of parameters to validate.
#' @param rules Named list mapping field name to the character vector of
#' allowed choices.
#' @param label Short human-readable label used in the error message
#' (e.g. `"MELD params"`).
#' @param hint Optional string appended to the error message. Defaults to
#' `NULL` (no hint).
#'
#' @return `TRUE` if all checked fields pass, otherwise a string of the form
#' `` "The element `<field>` in <label> is not one of the expected choices. <hint>" ``.
#'
#' @keywords internal
apply_choice_rules <- function(x, rules, label, hint = NULL) {
  res <- purrr::imap_lgl(x, \(val, name) {
    if (name %in% names(rules)) {
      checkmate::testChoice(val, rules[[name]])
    } else {
      TRUE
    }
  })
  if (all(res)) {
    return(TRUE)
  }
  broken <- names(res)[!res][1]
  msg <- sprintf(
    "The element `%s` in %s is not one of the expected choices.",
    broken,
    label
  )
  if (!is.null(hint)) {
    msg <- paste(msg, hint)
  }
  msg
}

## checks ----------------------------------------------------------------------

### others ---------------------------------------------------------------------

#' Check that files exist
#'
#' @description Checkmate extension for checking if files exist in the
#' directory.
#'
#' @param x String. The directory to check the files for.
#' @param file_names String. Vector of names of the expected files in this
#' directory.
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkFilesExist <- function(x, file_names) {
  res <- purrr::map(file_names, \(file) {
    checkmate::checkFileExists(file.path(x, file))
  })
  res <- purrr::keep(res, \(r) !is.logical(r))
  if (length(res) == 0) {
    return(TRUE)
  }
  res[[1]]
}

#' Assert that files exist
#'
#' @description Checkmate extension for asserting if files exist in the
#' directory.
#'
#' @inheritParams checkFilesExist
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertFileExists <- checkmate::makeAssertionFunction(checkFilesExist)

### correlation params ---------------------------------------------------------

#' Check correlation graph parameters
#'
#' @description Checkmate extension for checking the graph parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkCorGraphParams <- function(x) {
  res <- check_list_shape(
    x,
    c("epsilon", "min_cor", "fdr_threshold", "verbose")
  )
  if (!isTRUE(res)) {
    return(res)
  }

  apply_qtest_rules(
    x,
    list(
      epsilon = "R1",
      min_cor = "R1[0, 1]",
      fdr_threshold = "R1[0, 1]",
      verbose = "B1"
    ),
    label = "correlation graph params",
    hint = paste(
      "min_cor and fdr_threshold must be in [0, 1];",
      "epsilon must be a double; verbose must be a boolean."
    )
  )
}

#' Assert correlation graph parameters
#'
#' @description Checkmate extension for asserting the graph parameters.
#'
#' @inheritParams checkCorGraphParams
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertCorGraphParams <- checkmate::makeAssertionFunction(checkCorGraphParams)

### graph resolution -----------------------------------------------------------

#' Check resolution graph parameters
#'
#' @description Checkmate extension for checking the resolution parameters for
#' community detection with Leiden.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkGraphResParams <- function(x) {
  res <- check_list_shape(x, c("min_res", "max_res", "number_res"))
  if (!isTRUE(res)) {
    return(res)
  }

  apply_qtest_rules(
    x,
    list(min_res = "R1", max_res = "R1", number_res = "I1"),
    label = "resolution params",
    hint = "min_res and max_res must be doubles; number_res must be an integer."
  )
}

#' Assert resolution graph parameters
#'
#' @description Checkmate extension for asserting the resolution parameters for
#' community detection with Leiden.
#'
#' @inheritParams checkGraphResParams
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertGraphResParams <- checkmate::makeAssertionFunction(checkGraphResParams)

### ica ------------------------------------------------------------------------

#### ica params ----------------------------------------------------------------

#' Check ICA parameters
#'
#' @description Checkmate extension for checking the ICA parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkIcaParams <- function(x) {
  res <- check_list_shape(x, c("maxit", "alpha", "max_tol", "verbose"))
  if (!isTRUE(res)) {
    return(res)
  }

  apply_qtest_rules(
    x,
    list(
      maxit = "I1",
      alpha = "R1[1, 2]",
      max_tol = "R1(0, 1)",
      verbose = "B1"
    ),
    label = "ICA params",
    hint = paste(
      "maxit must be an integer; alpha must be in [1, 2];",
      "0 < max_tol < 1; verbose must be a boolean."
    )
  )
}

#' Assert ICA parameters
#'
#' @description Checkmate extension for asserting the ICA parameters
#'
#' @inheritParams checkIcaParams
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertIcaParams <- checkmate::makeAssertionFunction(checkIcaParams)

#### ica components ------------------------------------------------------------

#' Check ICA no of component parameters
#'
#' @description Checkmate extension for checking the ICA number of component
#' parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkIcaNcomps <- function(x) {
  res <- check_list_shape(x, c("max_no_comp", "steps", "custom_seq"))
  if (!isTRUE(res)) {
    return(res)
  }

  apply_qtest_rules(
    x,
    list(
      max_no_comp = "I1",
      steps = "I1",
      custom_seq = c("0", "I+")
    ),
    label = "ICA n-components params",
    hint = paste(
      "max_no_comp and steps must be integers;",
      "custom_seq must be NULL or a vector of integers."
    )
  )
}

#' Assert ICA no of component parameters
#'
#' @description Checkmate extension for asserting the ICA number of component
#' parameters.
#'
#' @inheritParams checkIcaNcomps
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertIcaNcomps <- checkmate::makeAssertionFunction(checkIcaNcomps)

#### ica randomisation ---------------------------------------------------------

#' Check ICA randomisation parameters
#'
#' @description Checkmate extension for checking the ICA randomisation
#' parameters for a version of stabilised ICA.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkIcaIterParams <- function(x) {
  res <- check_list_shape(x, c("cross_validate", "random_init", "folds"))
  if (!isTRUE(res)) {
    return(res)
  }

  apply_qtest_rules(
    x,
    list(
      cross_validate = "B1",
      random_init = "I1",
      folds = "I1"
    ),
    label = "ICA randomisation params",
    hint = paste(
      "random_init and folds must be integers;",
      "cross_validate must be a boolean."
    )
  )
}

#' Assert ICA randomisation parameters
#'
#' @description Checkmate extension for asserting the ICA randomisation
#' parameters for a version of stabilised ICA.
#'
#' @inheritParams checkIcaIterParams
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertIcaIterParams <- checkmate::makeAssertionFunction(checkIcaIterParams)

### community detections -------------------------------------------------------

#' Check community detection parameters
#'
#' @description Checkmate extension for checking the community detection
#' parameters for identifying genetically privileged communities.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkCommunityParams <- function(x) {
  res <- check_list_shape(
    x,
    c(
      "max_nodes",
      "min_nodes",
      "min_seed_nodes",
      "initial_res",
      "threshold_type",
      "network_threshold",
      "pval_threshold"
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  res <- apply_choice_rules(
    x,
    list(threshold_type = c("prop_based", "pval_based")),
    label = "community params"
  )
  if (!isTRUE(res)) {
    return(res)
  }

  apply_qtest_rules(
    x,
    list(
      max_nodes = sprintf("I1[%i,)", x$min_nodes),
      min_nodes = "I1",
      min_seed_nodes = "I1",
      initial_res = "N1",
      network_threshold = "N1(0, 1]",
      pval_threshold = "N1(0, 1]"
    ),
    label = "community params",
    hint = paste(
      "min_nodes, max_nodes and min_seed_nodes must be integers",
      "(with max_nodes >= min_nodes); initial_res must be a double;",
      "network_threshold and pval_threshold must be in (0, 1]."
    )
  )
}

#' Assert community detection parameter
#'
#' @description Checkmate extension for asserting the community detection
#' parameters for identifying genetically privileged communities.
#'
#' @inheritParams checkCommunityParams
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertCommunityParams <- checkmate::makeAssertionFunction(checkCommunityParams)

### gsea -----------------------------------------------------------------------

#' Check GSEA parameters
#'
#' @description Checkmate extension for checking the gene set enrichment
#' analysis (GSEA) parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkGSEAParams <- function(x) {
  res <- check_list_shape(
    x,
    c("min_size", "max_size", "gsea_param", "sample_size", "eps")
  )
  if (!isTRUE(res)) {
    return(res)
  }

  apply_qtest_rules(
    x,
    list(
      min_size = "I1[3,)",
      max_size = "I1[4,)",
      gsea_param = "N1",
      sample_size = "I1",
      eps = "N1"
    ),
    label = "GSEA params",
    hint = paste(
      "min_size and max_size must be integers (with max_size > min_size",
      "and min_size >= 3); gsea_param must be a double;",
      "sample_size must be an integer; eps must be a float."
    )
  )
}

#' Assert GSEA parameter
#'
#' @description Checkmate extension for asserting the gene set enrichment
#' analysis (GSEA) parameters.
#'
#' @inheritParams checkGSEAParams
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertGSEAParams <- checkmate::makeAssertionFunction(checkGSEAParams)

### gsva -----------------------------------------------------------------------

#' Check GSVA parameters
#'
#' @description Checkmate extension for checking the gene set variation analysis
#' (GSVA) parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkGSVAParams <- function(x) {
  res <- check_list_shape(
    x,
    c("tau", "min_size", "max_size", "max_diff", "abs_rank")
  )
  if (!isTRUE(res)) {
    return(res)
  }

  apply_qtest_rules(
    x,
    list(
      tau = "N1",
      min_size = "I1[3,)",
      max_size = "I1[4,)",
      max_diff = "B1",
      abs_rank = "B1"
    ),
    label = "GSVA params",
    hint = paste(
      "min_size and max_size must be integers (with max_size > min_size",
      "and min_size >= 3); tau must be a double;",
      "max_diff and abs_rank must be booleans."
    )
  )
}

#' Assert GSVA parameter
#'
#' @description Checkmate extension for asserting the gene set variation
#' analysis (GSVA) parameters.
#'
#' @inheritParams checkGSVAParams
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertGSVAParams <- checkmate::makeAssertionFunction(checkGSVAParams)

### ssgsea ---------------------------------------------------------------------

#' Check ssGSEA parameters
#'
#' @description Checkmate extension for checking single sample gene set
#' enrichment analysis parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkSingleSampleGSEAparams <- function(x) {
  res <- check_list_shape(
    x,
    c("alpha", "min_size", "max_size", "normalise")
  )
  if (!isTRUE(res)) {
    return(res)
  }

  apply_qtest_rules(
    x,
    list(
      alpha = "N1(0,1)",
      min_size = "I1[3,)",
      max_size = "I1[4,)",
      normalise = "B1"
    ),
    label = "ssGSEA params",
    hint = paste(
      "min_size and max_size must be integers (with max_size > min_size",
      "and min_size >= 3); alpha must be a double in (0, 1);",
      "normalise must be a boolean."
    )
  )
}

#' Assert ssGSEA parameter
#'
#' @description Checkmate extension for asserting single sample gene set
#' enrichment analysis parameters.
#'
#' @inheritParams checkSingleSampleGSEAparams
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertSingleSampleGSEAparams <- checkmate::makeAssertionFunction(
  checkSingleSampleGSEAparams
)

### coremo ---------------------------------------------------------------------

#' Check CoReMo parameters
#'
#' @description Checkmate extension for checking the CoReMo parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkCoReMoParams <- function(x) {
  res <- check_list_shape(
    x,
    c(
      "epsilon",
      "k_min",
      "k_max",
      "min_size",
      "junk_module_threshold",
      "rbf_func",
      "cor_method"
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  res <- apply_qtest_rules(
    x,
    list(
      epsilon = "N1",
      k_min = "I1",
      k_max = "I1",
      junk_module_threshold = "N1",
      min_size = c("I1", "0")
    ),
    label = "CoReMo params",
    hint = paste(
      "k_min and k_max must be integers; min_size must be an integer or NULL;",
      "junk_module_threshold and epsilon must be floats."
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  apply_choice_rules(
    x,
    list(
      rbf_func = c("gaussian", "inverse_quadratic", "bump"),
      cor_method = c("spearman", "pearson")
    ),
    label = "CoReMo params"
  )
}

#' Assert CoReMo parameter
#'
#' @description Checkmate extension for asserting the CoReMo parameters.
#'
#' @inheritParams checkCoReMoParams
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertCoReMoParams <- checkmate::makeAssertionFunction(checkCoReMoParams)

### dgrdl ----------------------------------------------------------------------

#' Check DGRDL parameters
#'
#' @description Checkmate extension for checking dual graph regularised
#' dictionary learning parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkDGRDLparams <- function(x) {
  res <- check_list_shape(
    x,
    c(
      "sparsity",
      "dict_size",
      "alpha",
      "beta",
      "max_iter",
      "k_neighbours",
      "admm_iter",
      "rho"
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  apply_qtest_rules(
    x,
    list(
      sparsity = "I1",
      dict_size = "I1",
      alpha = "N1",
      beta = "N1",
      max_iter = "I1",
      k_neighbours = "I1",
      admm_iter = "I1",
      rho = "N1"
    ),
    label = "DGRDL params",
    hint = paste(
      "sparsity, dict_size, max_iter, k_neighbours and admm_iter must be",
      "integers; alpha, beta and rho must be floats."
    )
  )
}

#' Assert DGRDL parameter
#'
#' @description Checkmate extension for asserting dual graph regularised
#' dictionary learning parameters.
#'
#' @inheritParams checkDGRDLparams
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertDGRDLparams <- checkmate::makeAssertionFunction(checkDGRDLparams)

### nmf (hals) -----------------------------------------------------------------

#' Check NMF HALS parameters
#'
#' @description Checkmate extension for checking the NMF HALS parameters for
#' single cell.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkNmfHals <- function(x) {
  res <- check_list_shape(
    x,
    c("max_iter", "tol", "eps", "check_every", "nmf_init")
  )
  if (!isTRUE(res)) {
    return(res)
  }

  res <- apply_qtest_rules(
    x,
    list(
      max_iter = "I1[1,)",
      tol = "N1(0,)",
      eps = "N1(0,)",
      check_every = "I1[1,)"
    ),
    label = "NMF HALS params",
    hint = paste(
      "max_iter and check_every must be positive integers;",
      "tol and eps must be positive numerics."
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  apply_choice_rules(
    x,
    list(nmf_init = c("nndsvd", "svd", "random")),
    label = "NMF HALS params"
  )
}

#' Assert NMF HALS parameters
#'
#' @description Checkmate extension for asserting the NMF HALS parameters for
#' single cell.
#'
#' @inheritParams checkNmfHals
#'
#' @param .var.name Name of the checked object to print in assertions.
#' @param add Collection to store assertion messages.
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertNmfHals <- checkmate::makeAssertionFunction(checkNmfHals)

### snf ------------------------------------------------------------------------

#' Check SNF parameters
#'
#' @description Checkmate extension for checking the SNF parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkSNFParams <- function(x) {
  res <- check_list_shape(
    x,
    c("k", "t", "mu", "alpha", "normalise", "distance_metric")
  )
  if (!isTRUE(res)) {
    return(res)
  }

  res <- apply_qtest_rules(
    x,
    list(
      k = "I1",
      t = "I1",
      mu = "N1[0,1]",
      alpha = "N1",
      normalise = "B1"
    ),
    label = "SNF params",
    hint = paste(
      "k and t must be positive integers; mu must be a float in [0, 1];",
      "alpha must be a float; normalise must be a boolean."
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  apply_choice_rules(
    x,
    list(distance_metric = c("euclidean", "manhattan", "canberra", "cosine")),
    label = "SNF params"
  )
}

#' Assert SNF parameter
#'
#' @description Checkmate extension for asserting the SNF parameters.
#'
#' @inheritParams checkSNFParams
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertSNFParams <- checkmate::makeAssertionFunction(checkSNFParams)

#' @rdname checkSNFParams
#'
#' @export
#'
#' @keywords internal
testSNFParams <- checkmate::makeTestFunction(checkSNFParams)

### cistarget ------------------------------------------------------------------

#' Check CisTarget parameters
#'
#' @description Checkmate extension for checking CisTarget parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkCistargetParams <- function(x) {
  res <- check_list_shape(
    x,
    c(
      "auc_threshold",
      "nes_threshold",
      "rcc_method",
      "high_conf_cats",
      "low_conf_cats"
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  res <- apply_qtest_rules(
    x,
    list(
      auc_threshold = "N1[0,1]",
      nes_threshold = "N1",
      rcc_method = "S1",
      high_conf_cats = "S+",
      low_conf_cats = "S+"
    ),
    label = "CisTarget params",
    hint = paste(
      "auc_threshold must be numeric [0, 1]; nes_threshold must be numeric;",
      "rcc_method must be a single string;",
      "high_conf_cats and low_conf_cats must be character vectors."
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  apply_choice_rules(
    x,
    list(rcc_method = c("approx", "icistarget")),
    label = "CisTarget params"
  )
}

#' Assert CisTarget parameters
#'
#' @description Checkmate extension for asserting the CisTarget parameters.
#'
#' @inheritParams checkCistargetParams
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertCistargetParams <- checkmate::makeAssertionFunction(checkCistargetParams)

### graph label propagation ----------------------------------------------------

#' Check label propagation parameters
#'
#' @description Checkmate extension for checking label propagation parameters.
#'
#' @param x The list to check/assert.
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkLabelPropParams <- function(x) {
  res <- check_list_shape(
    x,
    c("alpha", "iter", "tolerance", "symmetrise", "symmetry_strategy")
  )
  if (!isTRUE(res)) {
    return(res)
  }

  res <- apply_qtest_rules(
    x,
    list(
      alpha = "R1[0, 1]",
      iter = "I1[1,]",
      tolerance = "R1",
      symmetrise = "B1",
      symmetry_strategy = "S1"
    ),
    label = "label propagation params",
    hint = paste(
      "alpha must be in [0, 1]; iter must be a positive integer;",
      "tolerance must be a double; symmetrise must be a boolean;",
      "symmetry_strategy must be a string."
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  if (
    !is.null(x[["max_hops"]]) && !checkmate::qtest(x[["max_hops"]], "I1[0,]")
  ) {
    return("`max_hops` must be a positive integer or NULL.")
  }

  TRUE
}

#' Assert label propagation parameters
#'
#' @description Checkmate extension for asserting label propagation parameters.
#'
#' @inheritParams checkLabelPropParams
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertLabelPropParams <- checkmate::makeAssertionFunction(checkLabelPropParams)

### module membership ----------------------------------------------------------

#' Check module membership parameters
#'
#' @description Checkmate extension for checking the module membership
#' parameters, see [bixverse::params_module_membership()].
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkModuleMembershipParams <- function(x) {
  res <- check_list_shape(x, c("method", "cutoff", "fdr", "tails"))
  if (!isTRUE(res)) {
    return(res)
  }

  res <- apply_qtest_rules(
    x,
    list(cutoff = "N1(0,)", fdr = "N1(0,1]"),
    label = "module membership params",
    hint = "cutoff must be a positive float; fdr must sit in (0, 1]."
  )
  if (!isTRUE(res)) {
    return(res)
  }

  apply_choice_rules(
    x,
    list(
      method = c("zscore", "fdr"),
      tails = c("auto", "upper", "both")
    ),
    label = "module membership params"
  )
}

#' Assert module membership parameters
#'
#' @description Checkmate extension for asserting the module membership
#' parameters.
#'
#' @inheritParams checkModuleMembershipParams
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertModuleMembershipParams <- checkmate::makeAssertionFunction(
  checkModuleMembershipParams
)

### synthetic data (bulk) ------------------------------------------------------

#' Check synthetic bulk RNAseq parameters
#'
#' @description Checkmate extension for checking the synthetic bulk RNAseq
#' parameters, see [bixverse::params_synthetic_bulk_rnaseq()].
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkSyntheticBulkParams <- function(x) {
  res <- check_list_shape(
    x,
    c(
      "num_samples",
      "num_genes",
      "module_sizes",
      "generator",
      "seed",
      "mean_exp_gamma_shape",
      "mean_exp_gamma_scale",
      "disp_intercept",
      "disp_slope",
      "noise_std",
      "factor_std",
      "factor_shape",
      "factor_scale",
      "loading_mu",
      "loading_sigma",
      "hub_percentile"
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  res <- apply_qtest_rules(
    x,
    list(
      num_samples = "I1[1,)",
      num_genes = "I1[1,)",
      module_sizes = "I*",
      seed = "I1[0,)",
      mean_exp_gamma_shape = "N1(0,)",
      mean_exp_gamma_scale = "N1(0,)",
      disp_intercept = "N1(0,)",
      disp_slope = "N1(0,)",
      noise_std = "N1[0,)",
      factor_std = "N1(0,)",
      factor_shape = "N1(0,)",
      factor_scale = "N1(0,)",
      loading_mu = "N1",
      loading_sigma = "N1(0,)",
      hub_percentile = "N1(0,1]"
    ),
    label = "synthetic bulk params",
    hint = paste(
      "num_samples, num_genes and seed must be integers; module_sizes must be",
      "an integer vector (a double vector is silently ignored downstream);",
      "hub_percentile must sit in (0, 1]; the remaining distribution",
      "parameters must be positive floats."
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  res <- apply_choice_rules(
    x,
    list(
      generator = c(
        "hub_modular",
        "modular",
        "non_negative_factor",
        "non_gaussian_factor"
      )
    ),
    label = "synthetic bulk params"
  )
  if (!isTRUE(res)) {
    return(res)
  }

  if (sum(x[["module_sizes"]]) > x[["num_genes"]]) {
    return("The sum of `module_sizes` must not exceed `num_genes`.")
  }

  TRUE
}

#' Assert synthetic bulk RNAseq parameters
#'
#' @description Checkmate extension for asserting the synthetic bulk RNAseq
#' parameters.
#'
#' @inheritParams checkSyntheticBulkParams
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertSyntheticBulkParams <- checkmate::makeAssertionFunction(
  checkSyntheticBulkParams
)

#' Check bulk sparsification parameters
#'
#' @description Checkmate extension for checking the bulk sparsification
#' parameters, see [bixverse::params_bulk_sparsity()].
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkBulkSparsityParams <- function(x) {
  res <- check_list_shape(
    x,
    c("strategy", "target_library_size", "capture_efficiency_sigma", "seed")
  )
  if (!isTRUE(res)) {
    return(res)
  }

  res <- apply_qtest_rules(
    x,
    list(
      target_library_size = "N1(0,)",
      capture_efficiency_sigma = "N1(0,)",
      seed = "I1[0,)"
    ),
    label = "bulk sparsity params",
    hint = paste(
      "target_library_size and capture_efficiency_sigma must be positive",
      "floats; seed must be a positive integer."
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  apply_choice_rules(
    x,
    list(strategy = c("seq_depth")),
    label = "bulk sparsity params"
  )
}

#' Assert bulk sparsification parameters
#'
#' @description Checkmate extension for asserting the bulk sparsification
#' parameters.
#'
#' @inheritParams checkBulkSparsityParams
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertBulkSparsityParams <- checkmate::makeAssertionFunction(
  checkBulkSparsityParams
)

### single cell ----------------------------------------------------------------

#### general -------------------------------------------------------------------

#' Check kNN parameters
#'
#' @description Checkmate extension for checking kNN parameters.
#'
#' @param x The list to check/assert
#' @param required_params Character vector of required kNN parameter names
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkKnnParams <- function(x, required_params = NULL) {
  if (!is.null(required_params)) {
    res <- checkmate::checkNames(names(x), must.include = required_params)
    if (!isTRUE(res)) {
      return(res)
    }
  }

  res <- apply_qtest_rules(
    x,
    list(
      k = "I1[0,)",
      n_trees = "I1[1,)",
      search_budget = c("0", "I1[1,)"),
      m = "I1[1,)",
      ef_construction = "I1[1,)",
      ef_search = "I1[1,)",
      ef_budget = c("0", "I1[1,)"),
      n_list = c("0", "I1[1,)"),
      n_probe = c("0", "I1[1,)")
    ),
    label = "kNN params",
    hint = paste(
      "k must be >= 0;",
      "n_trees, m, ef_construction and ef_search must be >= 1;",
      "search_budget, ef_budget, n_list and n_probe must be NULL or >= 1."
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  res <- apply_qtest_rules(
    x,
    list(
      delta = "N1[0,1]",
      diversify_prob = "N1[0,1]"
    ),
    label = "kNN params",
    hint = "delta and diversify_prob must be in [0, 1]."
  )
  if (!isTRUE(res)) {
    return(res)
  }

  apply_choice_rules(
    x,
    list(
      knn_method = c(
        "annoy",
        "hnsw",
        "nndescent",
        "exhaustive",
        "ivf",
        "kmknn"
      ),
      ann_dist = c("euclidean", "cosine")
    ),
    label = "kNN params"
  )
}

#' Check FastCluster default parameters
#'
#' @description Checkmate extension for checking FastCluster parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkFastClusterDefaultParams <- function(x) {
  res <- checkmate::checkNames(
    names(x),
    must.include = c("km_type", "n_centroids", "kmeans_iters", "batch_size")
  )
  if (!isTRUE(res)) {
    return(res)
  }

  res <- apply_qtest_rules(
    x,
    list(
      batch_size = "I1",
      kmeans_iters = "I1",
      n_centroids = c("I1", "0")
    ),
    label = "FastCluster params",
    hint = paste(
      "batch_size and kmeans_iters must be integers;",
      "n_centroids must be an integer or NULL."
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  apply_choice_rules(
    x,
    list(km_type = c("minibatch", "standard")),
    label = "FastCluster params"
  )
}

#' Check k-means method parameters
#'
#' @description Checkmate extension for checking the run parameters of the
#' k-means clustering methods.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkKMeansParams <- function(x) {
  res <- checkmate::checkNames(
    names(x),
    must.include = c("k_means_iter", "k_means_init", "gemm", "hamerly")
  )
  if (!isTRUE(res)) {
    return(res)
  }

  res <- apply_qtest_rules(
    x,
    list(
      k_means_iter = "I1",
      gemm = c("0", "B1"),
      hamerly = c("0", "B1")
    ),
    label = "k-means params",
    hint = paste(
      "k_means_iter must be an integer;",
      "gemm and hamerly must be booleans or NULL."
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  apply_choice_rules(
    x,
    list(k_means_init = c("parallel", "random")),
    label = "k-means params"
  )
}

#### synthetic data ------------------------------------------------------------

##### rna ----------------------------------------------------------------------

#' Check synthetic data parameters
#'
#' @description Checkmate extension for checking the synthetic data
#' parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkScSyntheticData <- function(x) {
  res <- check_list_shape(
    x,
    c(
      "n_cells",
      "n_genes",
      "marker_genes",
      "n_batches",
      "batch_effect_strength",
      "n_samples",
      "sample_bias"
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  res <- apply_qtest_rules(
    x,
    list(
      n_cells = "I1",
      n_genes = "I1",
      n_batches = "I1",
      n_samples = c("0", "I1")
    ),
    label = "synthetic data params",
    hint = paste(
      "n_cells, n_genes and n_batches must be integers;",
      "n_samples must be an integer or NULL."
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  res <- checkmate::checkList(x$marker_genes, types = "list", names = "named")
  if (!isTRUE(res)) {
    return("marker_genes must be a named list of lists.")
  }

  res <- apply_choice_rules(
    x,
    list(batch_effect_strength = c("strong", "medium", "weak")),
    label = "synthetic data params"
  )
  if (!isTRUE(res)) {
    return(res)
  }

  if (!is.null(x[["sample_bias"]])) {
    res <- apply_choice_rules(
      x,
      list(sample_bias = c("even", "slightly_uneven", "very_uneven")),
      label = "synthetic data params"
    )
    if (!isTRUE(res)) {
      return(res)
    }
  }

  TRUE
}

#' Assert synthetic data parameters
#'
#' @description Checkmate extension for asserting the synthetic data
#' parameters.
#'
#' @inheritParams checkScSyntheticData
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertScSyntheticData <- checkmate::makeAssertionFunction(checkScSyntheticData)

##### adt ----------------------------------------------------------------------

#' Check synthetic ADT data parameters
#'
#' @description Checkmate extension for checking the synthetic ADT data
#' parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkScSyntheticDataAdt <- function(x) {
  res <- check_list_shape(
    x,
    c(
      "n_cells",
      "n_proteins",
      "marker_genes",
      "n_batches",
      "isotype_controls",
      "batch_effect_strength"
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  res <- apply_qtest_rules(
    x,
    list(
      n_cells = "I1",
      n_proteins = "I1",
      n_batches = "I1",
      isotype_controls = "I+"
    ),
    label = "synthetic ADT params",
    hint = paste(
      "n_cells, n_proteins and n_batches must be integers;",
      "isotype_controls must be a non-empty integer vector."
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  res <- checkmate::checkList(x$marker_genes, types = "list", names = "named")
  if (!isTRUE(res)) {
    return("marker_genes must be a named list of lists.")
  }

  apply_choice_rules(
    x,
    list(batch_effect_strength = c("strong", "medium", "weak")),
    label = "synthetic ADT params"
  )
}

#' Assert synthetic ADT data parameters
#'
#' @description Checkmate extension for asserting the synthetic ADT data
#' parameters.
#'
#' @inheritParams checkScSyntheticDataAdt
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertionFunction()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertScSyntheticDataAdt <- checkmate::makeAssertionFunction(
  checkScSyntheticDataAdt
)

#### io ------------------------------------------------------------------------

#' Check SC MTX load parameters
#'
#' @description Checkmate extension for checking MTX loading parameters
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkScMtxIO <- function(x) {
  res <- check_list_shape(
    x,
    c("path_mtx", "path_obs", "path_var", "cells_as_rows", "has_hdr")
  )
  if (!isTRUE(res)) {
    return(res)
  }

  files_ok <- purrr::map_lgl(c("path_mtx", "path_obs", "path_var"), \(n) {
    checkmate::testFileExists(x[[n]])
  })
  if (!all(files_ok)) {
    return(paste(
      "Some of the files specified in the config for mtx ingest do not exist.",
      "Please check the provided params."
    ))
  }

  apply_qtest_rules(
    x,
    list(cells_as_rows = "B1", has_hdr = "B1"),
    label = "MTX IO params",
    hint = "cells_as_rows and has_hdr must be booleans."
  )
}

#' Assert SC MTX load parameters
#'
#' @description Checkmate extension for asserting MTX loading parameters.
#'
#' @inheritParams checkScMtxIO
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertScMtxIO <- checkmate::makeAssertionFunction(checkScMtxIO)

#### qc ------------------------------------------------------------------------

#' Check SC minimum QC parameters
#'
#' @description Checkmate extension for checking the minimum QC parameters
#' in single cell.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkScMinQC <- function(x) {
  res <- check_list_shape(
    x,
    c("min_unique_genes", "min_lib_size", "min_cells", "target_size")
  )
  if (!isTRUE(res)) {
    return(res)
  }

  apply_qtest_rules(
    x,
    list(
      min_unique_genes = "I1",
      min_lib_size = "I1",
      min_cells = "I1",
      target_size = "N1"
    ),
    label = "single cell QC params",
    hint = paste(
      "min_unique_genes, min_lib_size and min_cells must be integers;",
      "target_size must be a float."
    )
  )
}

#' Assert SC minimum QC parameters
#'
#' @description Checkmate extension for asserting the minimum QC parameters
#' in single cell.
#'
#' @inheritParams checkScMinQC
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertScMinQC <- checkmate::makeAssertionFunction(checkScMinQC)

#### scrublet ------------------------------------------------------------------

#' Check Scrublet parameters
#'
#' @description Checkmate extension for checking the Scrublet parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkScScrublet <- function(x) {
  res <- check_list_shape(
    x,
    c(
      "log_transform",
      "mean_center",
      "normalise_variance",
      "target_size",
      "min_gene_var_pctl",
      "hvg_method",
      "loess_span",
      "clip_max",
      "n_bins",
      "binning_strategy",
      "sim_doublet_ratio",
      "expected_doublet_rate",
      "stdev_doublet_rate",
      "manual_threshold",
      "n_bins_histogram",
      "no_pcs",
      "random_svd"
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  res <- checkKnnParams(x[names(x) %in% KNN_PARAM_NAMES])
  if (!isTRUE(res)) {
    return(res)
  }

  res <- apply_qtest_rules(
    x,
    list(
      no_pcs = "I1[1,)",
      n_bins_histogram = "I1[10,)",
      n_bins = "I1[1,)",
      min_gene_var_pctl = "N1[0,1]",
      loess_span = "N1(0,)",
      sim_doublet_ratio = "N1(0,)",
      expected_doublet_rate = "N1[0,1]",
      stdev_doublet_rate = "N1[0,1]",
      target_size = "N1[0,)",
      log_transform = "B1",
      mean_center = "B1",
      normalise_variance = "B1",
      random_svd = "B1",
      clip_max = c("0", "N1(0,)"),
      manual_threshold = c("0", "N1[0,)")
    ),
    label = "Scrublet params",
    hint = paste(
      "no_pcs must be >= 1; n_bins_histogram must be >= 10; n_bins must be >= 1;",
      "min_gene_var_pctl, expected_doublet_rate and stdev_doublet_rate must be in [0, 1];",
      "loess_span and sim_doublet_ratio must be > 0; target_size must be >= 0;",
      "log_transform, mean_center, normalise_variance and random_svd must be booleans;",
      "clip_max and manual_threshold must be NULL or positive numerics."
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  apply_choice_rules(
    x,
    list(
      hvg_method = c("vst", "mvb", "dispersion"),
      binning_strategy = c("equal_width", "equal_frequency")
    ),
    label = "Scrublet params"
  )
}

#' Assert Scrublet parameters
#'
#' @description Checkmate extension for asserting the Scrublet parameters.
#'
#' @inheritParams checkScScrublet
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#'   to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#'   [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertScScrublet <- checkmate::makeAssertionFunction(checkScScrublet)

#### boost ---------------------------------------------------------------------

#' Check Boost parameters
#'
#' @description Checkmate extension for checking the Boost parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkScBoost <- function(x) {
  res <- check_list_shape(
    x,
    c(
      "log_transform",
      "mean_center",
      "normalise_variance",
      "target_size",
      "min_gene_var_pctl",
      "hvg_method",
      "loess_span",
      "clip_max",
      "n_bins",
      "binning_strategy",
      "boost_rate",
      "replace",
      "no_pcs",
      "random_svd",
      "resolution",
      "n_iters",
      "p_thresh",
      "voter_thresh",
      "fast_cluster"
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  res <- checkKnnParams(x[names(x) %in% KNN_PARAM_NAMES])
  if (!isTRUE(res)) {
    return(res)
  }

  res <- checkFastClusterDefaultParams(x)
  if (!isTRUE(res)) {
    return(res)
  }

  res <- apply_qtest_rules(
    x,
    list(
      no_pcs = "I1[1,)",
      n_iters = "I1[1,)",
      n_bins = "I1[1,)",
      min_gene_var_pctl = "N1[0,1]",
      loess_span = "N1(0,)",
      boost_rate = "N1[0,1]",
      resolution = "N1(0,)",
      p_thresh = "N1(0,)",
      voter_thresh = "N1[0,1]",
      target_size = "N1(0,)",
      log_transform = "B1",
      mean_center = "B1",
      normalise_variance = "B1",
      replace = "B1",
      random_svd = "B1",
      fast_cluster = "B1",
      clip_max = c("0", "N1(0,)")
    ),
    label = "Boost params",
    hint = paste(
      "no_pcs, n_bins and n_iters must be >= 1;",
      "min_gene_var_pctl, boost_rate and voter_thresh must be in [0, 1];",
      "loess_span, resolution and p_thresh must be > 0; target_size must be > 0;",
      "log_transform, mean_center, normalise_variance, replace, random_svd and",
      "fast_cluster must be booleans;",
      "clip_max must be NULL or a positive numeric."
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  apply_choice_rules(
    x,
    list(
      hvg_method = c("vst", "mvb", "dispersion"),
      binning_strategy = c("equal_width", "equal_frequency")
    ),
    label = "Boost params"
  )
}

#' Assert Boost parameters
#'
#' @description Checkmate extension for asserting the Boost parameters.
#'
#' @inheritParams checkScBoost
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#'   to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#'   [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertScBoost <- checkmate::makeAssertionFunction(checkScBoost)

#### scdblfinder ---------------------------------------------------------------

#' Check scDblFinder parameters
#'
#' @description Checkmate extension for checking scDblFinder parameters.
#'
#' @param x The list to check/assert.
#'
#' @return `TRUE` if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkScDblFinder <- function(x) {
  res <- check_list_shape(
    x,
    c(
      "log_transform",
      "mean_center",
      "normalise_variance",
      "target_size",
      "n_genes",
      "no_pcs",
      "random_svd",
      "doublet_ratio",
      "heterotypic_bias",
      "cluster_resolution",
      "cluster_iters",
      "fast_cluster",
      "n_iterations",
      "n_trees",
      "max_depth",
      "learning_rate",
      "min_samples_leaf",
      "subsample_rate",
      "cv_folds",
      "cv_early_stop",
      "se_fraction",
      "include_pcs",
      "cxds_genes"
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  res <- checkKnnParams(x[names(x) %in% KNN_PARAM_NAMES])
  if (!isTRUE(res)) {
    return(res)
  }

  res <- checkFastClusterDefaultParams(x)
  if (!isTRUE(res)) {
    return(res)
  }

  apply_qtest_rules(
    x,
    list(
      n_genes = "I1[1,)",
      no_pcs = "I1[1,)",
      cluster_iters = "I1[1,)",
      n_iterations = "I1[1,)",
      n_trees = "I1[1,)",
      max_depth = "I1[1,)",
      min_samples_leaf = "I1[1,)",
      cv_folds = "I1[2,)",
      cv_early_stop = "I1[1,)",
      include_pcs = "I1[1,)",
      doublet_ratio = "N1(0,)",
      heterotypic_bias = "N1[0,1]",
      cluster_resolution = "N1(0,)",
      learning_rate = "N1(0,)",
      subsample_rate = "N1(0,1]",
      se_fraction = "N1[0,)",
      target_size = "N1(0,)",
      log_transform = "B1",
      mean_center = "B1",
      normalise_variance = "B1",
      random_svd = "B1",
      fast_cluster = "B1",
      expected_doublet_rate = c("0", "N1(0,1]"),
      manual_threshold = c("0", "N1[0,)"),
      cxds_genes = c("0", "I1")
    ),
    label = "scDblFinder params",
    hint = "Please check the documentation for the expected type/range."
  )
}

#' Assert scDblFinder parameters
#'
#' @description Checkmate extension for asserting scDblFinder parameters.
#'
#' @inheritParams checkScDblFinder
#'
#' @param .var.name Name of the checked object to print in assertions.
#' @param add Collection to store assertion messages.
#'
#' @return Invisibly returns the checked object if successful.
#'
#' @keywords internal
assertScDblFinder <- checkmate::makeAssertionFunction(checkScDblFinder)

#### hvg -----------------------------------------------------------------------

#' Check HVG selection parameters
#'
#' @description Checkmate extension for checking the HVG parameters for single
#' cell.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkScHvg <- function(x) {
  res <- check_list_shape(
    x,
    c("method", "loess_span", "num_bin", "bin_method")
  )
  if (!isTRUE(res)) {
    return(res)
  }

  res <- apply_qtest_rules(
    x,
    list(
      loess_span = "N1[0.1, 1]",
      num_bin = "I1"
    ),
    label = "HVG params",
    hint = "loess_span must be in [0.1, 1]; num_bin must be an integer."
  )
  if (!isTRUE(res)) {
    return(res)
  }

  apply_choice_rules(
    x,
    list(
      method = c("vst", "meanvarbin", "dispersion"),
      bin_method = c("equal_width", "equal_freq")
    ),
    label = "HVG params"
  )
}

#' Assert HVG selection parameters
#'
#' @description Checkmate extension for checking the HVG parameters for single
#' cell.
#'
#' @inheritParams checkScHvg
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertScHvg <- checkmate::makeAssertionFunction(checkScHvg)

#### pca -----------------------------------------------------------------------

#' Check PCA parameters for single cell
#'
#' @description Checkmate extension for checking the PCA parameters for single
#' cell.
#'
#' @param x The list to check/assert.
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkScPca <- function(x) {
  res <- check_list_shape(
    x,
    c("mean_center", "normalise_variance", "randomised", "clr", "size_factor")
  )
  if (!isTRUE(res)) {
    return(res)
  }

  apply_qtest_rules(
    x,
    list(
      mean_center = "B1",
      normalise_variance = "B1",
      randomised = "B1",
      clr = "B1",
      size_factor = "N1"
    ),
    label = "single cell PCA params",
    hint = paste(
      "mean_center, normalise_variance, randomised and clr must be single",
      "booleans; size_factor must be a single numeric."
    )
  )
}

#' Assert PCA parameters for single cell
#'
#' @description Checkmate extension for asserting the PCA parameters for single
#' cell.
#'
#' @inheritParams checkScPca
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertScPca <- checkmate::makeAssertionFunction(checkScPca)

#### knn -----------------------------------------------------------------------

#' Check single cell kNN parameters
#'
#' @description Checkmate extension for checking the single cell kNN
#' parameters.
#'
#' @param x The list to check/assert.
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkScKnn <- function(x) {
  res <- check_list_shape(
    x,
    c(
      "k",
      "knn_method",
      "ann_dist",
      "n_trees",
      "search_budget",
      "delta",
      "diversify_prob",
      "ef_budget",
      "m",
      "ef_construction",
      "ef_search",
      "n_list",
      "n_probe"
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  res <- apply_qtest_rules(
    x,
    list(
      k = "I1[1,)",
      n_trees = "I1[1,)",
      delta = "N1(0,)",
      diversify_prob = "N1[0,1]",
      m = "I1[1,)",
      ef_construction = "I1[1,)",
      ef_search = "I1[1,)",
      search_budget = c("0", "I1[1,)"),
      ef_budget = c("0", "I1[1,)"),
      n_list = c("0", "I1[1,)"),
      n_probe = c("0", "I1[1,)")
    ),
    label = "kNN params",
    hint = paste(
      "k, n_trees, m, ef_construction and ef_search must be positive integers;",
      "delta must be a positive numeric;",
      "diversify_prob must be a numeric in [0, 1];",
      "search_budget, ef_budget, n_list and n_probe must be NULL or positive integers."
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  apply_choice_rules(
    x,
    list(
      knn_method = c(
        "kmknn",
        "hnsw",
        "annoy",
        "nndescent",
        "ivf",
        "exhaustive"
      ),
      ann_dist = c("cosine", "euclidean")
    ),
    label = "kNN params"
  )
}

#' Assert single cell kNN parameters
#'
#' @description Checkmate extension for asserting the single cell kNN
#' parameters.
#'
#' @inheritParams checkScKnn
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertionFunction()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertScKnn <- checkmate::makeAssertionFunction(checkScKnn)

#### neighbours ----------------------------------------------------------------

#' Check neighbour generation parameters
#'
#' @description Checkmate extension for checking the neighbour generation
#' parameters for single cell.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkScNeighbours <- function(x) {
  res <- check_list_shape(
    x,
    c("full_snn", "pruning", "snn_similarity", "ann_dist")
  )
  if (!isTRUE(res)) {
    return(res)
  }

  res <- checkKnnParams(x[names(x) %in% KNN_PARAM_NAMES])
  if (!isTRUE(res)) {
    return(res)
  }

  res <- apply_qtest_rules(
    x,
    list(full_snn = "B1", pruning = "N1[0, 1]"),
    label = "neighbour params",
    hint = "full_snn must be a boolean; pruning must be in [0, 1]."
  )
  if (!isTRUE(res)) {
    return(res)
  }

  apply_choice_rules(
    x,
    list(snn_similarity = c("rank", "jaccard")),
    label = "neighbour params"
  )
}

#' Assert neighbour generation parameters
#'
#' @description Checkmate extension for assert the neighbour generation
#' parameters for single cell.
#'
#' @inheritParams checkScNeighbours
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertScNeighbours <- checkmate::makeAssertionFunction(checkScNeighbours)

#### cells in object -----------------------------------------------------------

#' Check that the cell name exists in the object
#'
#' @description Checkmate extension for checking if the prodivided cell names
#' exist in the object.
#'
#' @param x The `SingleCells` or `SingleCellsSubset` object to check/assert.
#' @param cell_names String. The provided cell names.
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkCellsExist <- function(x, cell_names) {
  res <- checkmate::checkMultiClass(
    x,
    c("bixverse::SingleCells", "bixverse::SingleCellsSubset")
  )
  if (!isTRUE(res)) {
    return(res)
  }
  res <- checkmate::qtest(cell_names, "S+")
  if (!isTRUE(res)) {
    return("The cell names need to be a string vector.")
  }
  if (!all(cell_names %in% get_cell_names(x))) {
    return(paste(
      "Some of the provided cell names do not exist in the object.",
      "Please check."
    ))
  }
  TRUE
}

#' Assert neighbour generation parameters
#'
#' @description Checkmate extension for asserting if the prodivided cell names
#  exist in the object.
#'
#' @inheritParams checkCellsExist
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertCellsExist <- checkmate::makeAssertionFunction(checkCellsExist)

#### meta cells ----------------------------------------------------------------

#' Check (bootstrapped) metacell generation parameters
#'
#' @description Checkmate extension for checking the (bootstrapped) metacell
#' generation parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkScBootstrappedMetacells <- function(x) {
  res <- check_list_shape(
    x,
    c("max_shared", "target_no_metacells", "max_iter")
  )
  if (!isTRUE(res)) {
    return(res)
  }

  res <- checkKnnParams(x[names(x) %in% KNN_PARAM_NAMES])
  if (!isTRUE(res)) {
    return(res)
  }

  apply_qtest_rules(
    x,
    list(
      max_shared = "I1[1,)",
      target_no_metacells = "I1[1,)",
      max_iter = "I1[1,)"
    ),
    label = "bootstrapped metacell params",
    hint = "max_shared, target_no_metacells and max_iter must be integers >= 1."
  )
}

#' Assert (bootstrapped) metacell generation parameters
#'
#' @description Checkmate extension for assert the metacell generation
#' parameters.
#'
#' @inheritParams checkScBootstrappedMetacells
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertScBootstrappedMetacells <- checkmate::makeAssertionFunction(
  checkScBootstrappedMetacells
)

#### seacells ------------------------------------------------------------------

#' Check SEACells parameters
#'
#' @description Checkmate extension for checking the SEACells parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkScSeacells <- function(x) {
  res <- check_list_shape(
    x,
    c(
      "n_sea_cells",
      "max_fw_iters",
      "convergence_epsilon",
      "max_iter",
      "min_iter",
      "greedy_threshold",
      "graph_building",
      "pruning",
      "pruning_threshold",
      "n_landmarks"
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  res <- checkKnnParams(x[names(x) %in% KNN_PARAM_NAMES])
  if (!isTRUE(res)) {
    return(res)
  }

  apply_qtest_rules(
    x,
    list(
      n_sea_cells = "I1[1,)",
      max_fw_iters = "I1[1,)",
      max_iter = "I1[1,)",
      min_iter = "I1[1,)",
      greedy_threshold = "I1[1,)",
      convergence_epsilon = "N1",
      pruning_threshold = "N1",
      pruning = "B1",
      graph_building = "S1",
      n_landmarks = c("0", "I1")
    ),
    label = "SEACells params",
    hint = paste(
      "n_sea_cells, max_fw_iters, max_iter, min_iter and greedy_threshold",
      "must be integers >= 1;",
      "convergence_epsilon and pruning_threshold must be numeric;",
      "pruning must be a boolean; graph_building must be a string;",
      "n_landmarks must be an integer or NULL."
    )
  )
}

#' Assert SEACells parameters
#'
#' @description Checkmate extension for asserting the SEACells parameters.
#'
#' @inheritParams checkScSeacells
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertScSeacells <- checkmate::makeAssertionFunction(checkScSeacells)

#### supercells ----------------------------------------------------------------

#' Check SuperCell parameters
#'
#' @description Checkmate extension for checking the SuperCell parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkScSupercell <- function(x) {
  res <- check_list_shape(
    x,
    c("walk_length", "graining_factor", "use_kernel", "k_ith", "max_support")
  )
  if (!isTRUE(res)) {
    return(res)
  }

  res <- checkKnnParams(x[names(x) %in% KNN_PARAM_NAMES])
  if (!isTRUE(res)) {
    return(res)
  }

  apply_qtest_rules(
    x,
    list(
      walk_length = "I1[1,)",
      k_ith = c("I1", "0"),
      max_support = c("I1[1,)", "0"),
      graining_factor = "N1",
      use_kernel = "B1"
    ),
    label = "SuperCell params",
    hint = paste(
      "walk_length must be an integer >= 1;",
      "k_ith and max_support must be an integer or NULL;",
      "graining_factor must be numeric; use_kernel must be a boolean."
    )
  )
}

#' Assert SuperCell parameters
#'
#' @description Checkmate extension for asserting the SuperCell parameters.
#'
#' @inheritParams checkScSupercell
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertScSupercell <- checkmate::makeAssertionFunction(checkScSupercell)

#### bbknn ---------------------------------------------------------------------

#' Check BBKNN parameters
#'
#' @description Checkmate extension for checking the BBKNN parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkScBbknn <- function(x) {
  res <- check_list_shape(
    x,
    c(
      "neighbours_within_batch",
      "set_op_mix_ratio",
      "local_connectivity",
      "trim"
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  res <- checkKnnParams(x[names(x) %in% KNN_PARAM_NAMES])
  if (!isTRUE(res)) {
    return(res)
  }

  apply_qtest_rules(
    x,
    list(
      neighbours_within_batch = "I1[1,)",
      trim = c("0", "I1[1,)"),
      set_op_mix_ratio = "N1[0,1]",
      local_connectivity = "N1"
    ),
    label = "BBKNN params",
    hint = paste(
      "neighbours_within_batch must be an integer >= 1;",
      "trim must be NULL or an integer >= 1;",
      "set_op_mix_ratio must be a numeric in [0, 1];",
      "local_connectivity must be numeric."
    )
  )
}

#' Assert BBKNN parameters
#'
#' @description Checkmate extension for asserting the BBKNN parameters.
#'
#' @inheritParams checkScBbknn
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertScBbknn <- checkmate::makeAssertionFunction(checkScBbknn)

#### fastMNN -------------------------------------------------------------------

#' Check fastMNN parameters
#'
#' @description Checkmate extension for checking the fastMNN parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkScFastmnn <- function(x) {
  res <- check_list_shape(
    x,
    c(
      "ndist",
      "cos_norm",
      "no_pcs",
      "randomised",
      "sparse_svd",
      "mean_center",
      "normalise_variance",
      "clr",
      "size_factor"
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  res <- checkKnnParams(x[names(x) %in% KNN_PARAM_NAMES])
  if (!isTRUE(res)) {
    return(res)
  }

  apply_qtest_rules(
    x,
    list(
      no_pcs = "I1[1,)",
      ndist = "N1(0,)",
      size_factor = "N1",
      cos_norm = "B1",
      randomised = "B1",
      sparse_svd = "B1",
      mean_center = "B1",
      normalise_variance = "B1",
      clr = "B1"
    ),
    label = "fastMNN params",
    hint = paste(
      "no_pcs must be an integer >= 1; ndist must be a positive numeric;",
      "size_factor must be numeric;",
      "cos_norm, randomised, sparse_svd, mean_center, normalise_variance",
      "and clr must be booleans."
    )
  )
}

#' Assert fastMNN parameters
#'
#' @description Checkmate extension for asserting the fastMNN parameters.
#'
#' @inheritParams checkScFastmnn
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertScFastmnn <- checkmate::makeAssertionFunction(checkScFastmnn)

#### seurat CCA ----------------------------------------------------------------

#' Check Seurat CCA parameters
#'
#' @description Checkmate extension for checking the Seurat CCA parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkScSeuratCca <- function(x) {
  res <- checkmate::checkList(x)
  if (!isTRUE(res)) {
    return(res)
  }
  res <- checkmate::checkNames(
    names(x),
    must.include = c(
      "num_cc",
      "dims",
      "k_anchor",
      "k_filter",
      "k_score",
      "k_weight",
      "n_top_features",
      "l2_norm",
      "sd",
      "randomised",
      "mean_center",
      "normalise_variance",
      "clr",
      "size_factor"
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }
  # knn
  knn_params <- x[names(x) %in% KNN_PARAM_NAMES]
  res <- checkKnnParams(knn_params)
  if (!isTRUE(res)) {
    return(res)
  }
  # Check integer parameters. k_filter may be 0 to disable the filter step.
  int_params <- c(
    "num_cc",
    "dims",
    "k_anchor",
    "k_score",
    "k_weight",
    "n_top_features"
  )
  res <- purrr::map_lgl(int_params, \(param) {
    checkmate::qtest(x[[param]], "I1[1,)")
  })
  if (!isTRUE(all(res))) {
    broken_param <- int_params[which(!res)][1]
    return(sprintf("%s needs to be an integer >= 1.", broken_param))
  }
  res <- checkmate::qtest(x[["k_filter"]], "I1[0,)")
  if (!isTRUE(res)) {
    return("k_filter needs to be an integer >= 0.")
  }
  # Check numeric parameters
  res <- checkmate::qtest(x[["sd"]], "N1(0,)")
  if (!isTRUE(res)) {
    return("sd needs to be a positive numeric.")
  }
  res <- checkmate::checkNumber(x[["size_factor"]])
  if (!isTRUE(res)) {
    return(res)
  }
  # Check logical parameters
  logical_params <- c(
    "l2_norm",
    "randomised",
    "mean_center",
    "normalise_variance",
    "clr"
  )
  res <- purrr::map_lgl(logical_params, \(param) {
    checkmate::qtest(x[[param]], "B1")
  })
  if (!isTRUE(all(res))) {
    broken_param <- logical_params[which(!res)][1]
    return(sprintf("%s needs to be logical.", broken_param))
  }

  return(TRUE)
}

#' Assert Seurat CCA parameters
#'
#' @description Checkmate extension for asserting the Seurat CCA parameters.
#'
#' @inheritParams checkScSeuratCca
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertScSeuratCca <- checkmate::makeAssertionFunction(checkScSeuratCca)

#### seurat rPCA ---------------------------------------------------------------

#' Check Seurat rPCA parameters
#'
#' @description Checkmate extension for checking the Seurat rPCA parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkScSeuratRpca <- function(x) {
  res <- checkmate::checkList(x)
  if (!isTRUE(res)) {
    return(res)
  }
  res <- checkmate::checkNames(
    names(x),
    must.include = c(
      "dims",
      "k_anchor",
      "k_score",
      "k_weight",
      "l2_norm",
      "sd",
      "randomised",
      "mean_center",
      "normalise_variance",
      "clr",
      "size_factor"
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }
  # knn
  knn_params <- x[names(x) %in% KNN_PARAM_NAMES]
  res <- checkKnnParams(knn_params)
  if (!isTRUE(res)) {
    return(res)
  }
  # Check integer parameters
  int_params <- c("dims", "k_anchor", "k_score", "k_weight")
  res <- purrr::map_lgl(int_params, \(param) {
    checkmate::qtest(x[[param]], "I1[1,)")
  })
  if (!isTRUE(all(res))) {
    broken_param <- int_params[which(!res)][1]
    return(sprintf("%s needs to be an integer >= 1.", broken_param))
  }
  # Check numeric parameters
  res <- checkmate::qtest(x[["sd"]], "N1(0,)")
  if (!isTRUE(res)) {
    return("sd needs to be a positive numeric.")
  }
  res <- checkmate::checkNumber(x[["size_factor"]])
  if (!isTRUE(res)) {
    return(res)
  }
  # Check logical parameters
  logical_params <- c(
    "l2_norm",
    "randomised",
    "mean_center",
    "normalise_variance",
    "clr"
  )
  res <- purrr::map_lgl(logical_params, \(param) {
    checkmate::qtest(x[[param]], "B1")
  })
  if (!isTRUE(all(res))) {
    broken_param <- logical_params[which(!res)][1]
    return(sprintf("%s needs to be logical.", broken_param))
  }

  return(TRUE)
}

#' Assert Seurat rPCA parameters
#'
#' @description Checkmate extension for asserting the Seurat rPCA parameters.
#'
#' @inheritParams checkScSeuratRpca
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertScSeuratRpca <- checkmate::makeAssertionFunction(checkScSeuratRpca)

#### VISION --------------------------------------------------------------------

#' Check VISION parameters
#'
#' @description Checkmate extension for checking the VISION parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkScVision <- function(x) {
  res <- check_list_shape(x, c("n_perm", "n_cluster"))
  if (!isTRUE(res)) {
    return(res)
  }

  res <- checkKnnParams(x[names(x) %in% KNN_PARAM_NAMES])
  if (!isTRUE(res)) {
    return(res)
  }

  apply_qtest_rules(
    x,
    list(
      n_perm = "I1[1,)",
      n_cluster = "I1[1,)"
    ),
    label = "VISION params",
    hint = "n_perm and n_cluster must be integers >= 1."
  )
}

#' Assert VISION parameters
#'
#' @description Checkmate extension for asserting the VISION parameters.
#'
#' @inheritParams checkScVision
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertScVision <- checkmate::makeAssertionFunction(checkScVision)

#### AUCell --------------------------------------------------------------------

#' Check AUCell parameters
#'
#' @description Checkmate extension for checking the AUCell parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkScAucell <- function(x) {
  res <- check_list_shape(x, c("auc_type", "max_rank", "standardise"))
  if (!isTRUE(res)) {
    return(res)
  }

  res <- apply_qtest_rules(
    x,
    list(
      max_rank = c("N1[1,)", "0"),
      standardise = "B1"
    ),
    label = "AUCell params",
    hint = paste(
      "max_rank must be NULL or a single numeric >= 1;",
      "standardise must be a single logical."
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  apply_choice_rules(
    x,
    list(auc_type = c("wilcox", "recovery", "ap")),
    label = "AUCell params"
  )
}

#' Assert AUCell parameters
#'
#' @description Checkmate extension for asserting the AUCell parameters.
#'
#' @inheritParams checkScAucell
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertScAucell <- checkmate::makeAssertionFunction(checkScAucell)

#### HotSpot -------------------------------------------------------------------

#' Check HotSpot parameters
#'
#' @description Checkmate extension for checking the HotSpot parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkScHotspot <- function(x) {
  res <- check_list_shape(x, c("model", "normalise"))
  if (!isTRUE(res)) {
    return(res)
  }

  res <- checkKnnParams(x[names(x) %in% KNN_PARAM_NAMES])
  if (!isTRUE(res)) {
    return(res)
  }

  res <- apply_qtest_rules(
    x,
    list(normalise = "B1"),
    label = "HotSpot params",
    hint = "normalise must be a boolean."
  )
  if (!isTRUE(res)) {
    return(res)
  }

  apply_choice_rules(
    x,
    list(model = c("danb", "bernoulli", "normal")),
    label = "HotSpot params"
  )
}

#' Assert HotSpot parameters
#'
#' @description Checkmate extension for asserting the HotSpot parameters.
#'
#' @inheritParams checkScHotspot
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertScHotspot <- checkmate::makeAssertionFunction(checkScHotspot)

#### miloR ---------------------------------------------------------------------

#' Check MiloR parameters
#'
#' @description Checkmate extension for checking the MiloR parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkScMiloR <- function(x) {
  res <- check_list_shape(
    x,
    c("prop", "k_refine", "refinement_strategy", "index_type")
  )
  if (!isTRUE(res)) {
    return(res)
  }

  res <- checkKnnParams(x[names(x) %in% KNN_PARAM_NAMES])
  if (!isTRUE(res)) {
    return(res)
  }

  res <- apply_qtest_rules(
    x,
    list(
      k_refine = "I1[1,)",
      prop = "N1(0,1)"
    ),
    label = "MiloR params",
    hint = "k_refine must be an integer >= 1; prop must be in (0, 1)."
  )
  if (!isTRUE(res)) {
    return(res)
  }

  apply_choice_rules(
    x,
    list(
      refinement_strategy = c("approximate", "bruteforce", "index"),
      index_type = c("nndescent", "ivf", "hnsw", "annoy")
    ),
    label = "MiloR params"
  )
}

#' Assert MiloR parameters
#'
#' @description Checkmate extension for asserting the MiloR parameters.
#'
#' @inheritParams checkScMiloR
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertScMiloR <- checkmate::makeAssertionFunction(checkScMiloR)

#### Harmony -------------------------------------------------------------------

#' Check Harmony parameters
#'
#' @description Checkmate extension for checking Harmony parameters.
#'
#' @param x The list to check/assert.
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkScHarmonyParams <- function(x) {
  res <- checkmate::checkList(x)
  if (!isTRUE(res)) {
    return(res)
  }

  res <- checkKMeansParams(x)
  if (!isTRUE(res)) {
    return(res)
  }

  res <- checkmate::checkNames(
    names(x),
    must.include = c(
      "k",
      "sigma",
      "theta",
      "lambda",
      "block_size",
      "max_iter_kmeans",
      "max_iter_harmony",
      "epsilon_kmeans",
      "epsilon_harmony",
      "window_size"
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  apply_qtest_rules(
    x,
    list(
      k = c("I1[1,)", "0"),
      max_iter_kmeans = "I1[1,)",
      max_iter_harmony = "I1[1,)",
      window_size = "I1[1,)",
      sigma = "N+[0,)",
      theta = "N+[0,)",
      lambda = "N+[0,)",
      block_size = "N1(0,1]",
      epsilon_kmeans = "N1(0,)",
      epsilon_harmony = "N1(0,)"
    ),
    label = "Harmony params",
    hint = paste(
      "max_iter_kmeans, max_iter_harmony and window_size must be integers >= 1;",
      "k must be NULL or an integer;",
      "sigma, theta and lambda must be numeric vectors with non-negative values;",
      "block_size must be in (0, 1];",
      "epsilon_kmeans and epsilon_harmony must be > 0."
    )
  )
}

#' Assert Harmony parameters
#'
#' @description Checkmate extension for asserting the Harmony parameters.
#'
#' @inheritParams checkScHarmonyParams
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertScHarmonyParams <- checkmate::makeAssertionFunction(checkScHarmonyParams)

#### Harmony (version 2) -------------------------------------------------------

#' Check Harmony v2 parameters
#'
#' @description Checkmate extension for checking Harmony v2 parameters.
#'
#' @param x The list to check/assert.
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkScHarmonyParamsV2 <- function(x) {
  res <- checkmate::checkList(x)
  if (!isTRUE(res)) {
    return(res)
  }

  res <- checkKMeansParams(x)
  if (!isTRUE(res)) {
    return(res)
  }

  res <- checkmate::checkNames(
    names(x),
    must.include = c(
      "k",
      "sigma",
      "theta",
      "lambda",
      "block_size",
      "max_iter_kmeans",
      "max_iter_harmony",
      "epsilon_kmeans",
      "epsilon_harmony",
      "window_size",
      "alpha",
      "tau",
      "batch_proportion_cutoff",
      "use_dynamic_lambda"
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  apply_qtest_rules(
    x,
    list(
      k = c("I1[1,)", "0"),
      max_iter_kmeans = "I1[1,)",
      max_iter_harmony = "I1[1,)",
      window_size = "I1[1,)",
      sigma = "N+[0,)",
      theta = "N+[0,)",
      lambda = "N+[0,)",
      block_size = "N1(0,1]",
      epsilon_kmeans = "N1(0,)",
      epsilon_harmony = "N1(0,)",
      alpha = "N1(0,1)",
      tau = "N1[0,)",
      batch_proportion_cutoff = "N1(0,)",
      use_dynamic_lambda = "B1"
    ),
    label = "Harmony v2 params",
    hint = paste(
      "max_iter_kmeans, max_iter_harmony and window_size must be integers >= 1;",
      "k must be NULL or an integer;",
      "sigma, theta and lambda must be numeric vectors with non-negative values;",
      "block_size must be in (0, 1];",
      "epsilon_kmeans, epsilon_harmony and batch_proportion_cutoff must be > 0;",
      "alpha must be in (0, 1); tau must be >= 0;",
      "use_dynamic_lambda must be a single logical."
    )
  )
}

#' Assert Harmony v2 parameters
#'
#' @description Checkmate extension for asserting the Harmony v2 parameters.
#'
#' @inheritParams checkScHarmonyParamsV2
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertScHarmonyParamsV2 <- checkmate::makeAssertionFunction(
  checkScHarmonyParamsV2
)

#### symphony ------------------------------------------------------------------

#' Check Symphony map parameters
#'
#' @description Checkmate extension for checking the Symphony query mapping
#' parameters.
#'
#' @param x The list to check/assert.
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkSymphonyMap <- function(x) {
  res <- check_list_shape(x, c("sigma", "lambda"))
  if (!isTRUE(res)) {
    return(res)
  }

  apply_qtest_rules(
    x,
    list(sigma = "N1[0,)", lambda = "N1[0,)"),
    label = "Symphony map params",
    hint = "sigma and lambda must be non-negative numerics."
  )
}

#' Assert Symphony map parameters
#'
#' @description Checkmate extension for asserting the Symphony query mapping
#' parameters.
#'
#' @inheritParams checkSymphonyMap
#'
#' @param .var.name Name of the checked object to print in assertions.
#' @param add Collection to store assertion messages.
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertSymphonyMap <- checkmate::makeAssertionFunction(checkSymphonyMap)

#### scenic --------------------------------------------------------------------

#' Check SCENIC parameters
#'
#' @description Checkmate extension for checking SCENIC parameters.
#'
#' @param x The list to check/assert.
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkScenicParams <- function(x) {
  res <- check_list_shape(
    x,
    c(
      "min_counts",
      "min_cells",
      "learner_type",
      "gene_batch_strategy",
      "n_pcs",
      "n_subsample",
      "min_samples_leaf",
      "n_features_split",
      "max_depth"
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  res <- apply_choice_rules(
    x,
    list(
      learner_type = c("randomforest", "extratrees", "grnboost2"),
      gene_batch_strategy = c("random", "correlated")
    ),
    label = "SCENIC params"
  )
  if (!isTRUE(res)) {
    return(res)
  }

  res <- apply_qtest_rules(
    x,
    list(
      min_counts = "I1[1,)",
      n_pcs = "I1[1,)",
      n_subsample = "I1[1,)",
      min_samples_leaf = "I1[1,)",
      n_features_split = "I1[0,)",
      max_depth = "I1[1,)",
      min_cells = "N1(0,1]"
    ),
    label = "SCENIC params",
    hint = paste(
      "min_counts, n_pcs, n_subsample, min_samples_leaf and max_depth must be",
      "integers >= 1; n_features_split must be an integer >= 0;",
      "min_cells must be in (0, 1]."
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  if (!is.null(x$gene_batch_size)) {
    if (!checkmate::qtest(x$gene_batch_size, "I1[1,)")) {
      return("gene_batch_size must be a positive integer or NULL.")
    }
  }

  # Learner-specific validation
  if (x$learner_type == "randomforest") {
    if (is.null(x$n_trees) || !checkmate::qtest(x$n_trees, "I1[1,)")) {
      return("n_trees must be a positive integer for randomforest.")
    }
    if (
      is.null(x$subsample_rate) ||
        !checkmate::qtest(x$subsample_rate, "N1(0,1]")
    ) {
      return("subsample_rate must be a numeric in (0, 1] for randomforest.")
    }
    if (is.null(x$bootstrap) || !checkmate::qtest(x$bootstrap, "B1")) {
      return("bootstrap must be a single logical for randomforest.")
    }
    if (
      !is.null(x$subsample_frac) &&
        !checkmate::qtest(x$subsample_frac, "N1(0,1]")
    ) {
      return("subsample_frac must be a numeric in (0, 1] or NULL.")
    }
  }

  if (x$learner_type == "extratrees") {
    if (is.null(x$n_trees) || !checkmate::qtest(x$n_trees, "I1[1,)")) {
      return("n_trees must be a positive integer for extratrees.")
    }
    if (
      is.null(x$n_thresholds) ||
        !checkmate::qtest(x$n_thresholds, "I1[1,)")
    ) {
      return("n_thresholds must be a positive integer for extratrees.")
    }
    if (
      !is.null(x$subsample_frac) &&
        !checkmate::qtest(x$subsample_frac, "N1(0,1]")
    ) {
      return("subsample_frac must be a numeric in (0, 1] or NULL.")
    }
  }

  if (x$learner_type == "grnboost2") {
    if (is.null(x$n_trees_max) || !checkmate::qtest(x$n_trees_max, "I1[1,)")) {
      return("n_trees_max must be a positive integer for grnboost2.")
    }
    if (
      is.null(x$learning_rate) ||
        !checkmate::qtest(x$learning_rate, "N1(0,1]")
    ) {
      return("learning_rate must be a numeric in (0, 1] for grnboost2.")
    }
    if (
      is.null(x$early_stop_window) ||
        !checkmate::qtest(x$early_stop_window, "I1[1,)")
    ) {
      return("early_stop_window must be a positive integer for grnboost2.")
    }
    if (
      is.null(x$subsample_rate) ||
        !checkmate::qtest(x$subsample_rate, "N1(0,1]")
    ) {
      return("subsample_rate must be a numeric in (0, 1] for grnboost2.")
    }
  }

  TRUE
}

#' Assert SCENIC parameters
#'
#' @description Checkmate extension for asserting SCENIC parameters.
#'
#' @inheritParams checkScenicParams
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertScenicParams <- checkmate::makeAssertionFunction(checkScenicParams)

#### fast clustering -----------------------------------------------------------

#' Check singe cell fast clustering parameters
#'
#' @description Checkmate extension for checking the fast clustering parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkScFastCluster <- function(x) {
  res <- check_list_shape(
    x,
    c(
      "kmeans_iters",
      "batch_size",
      "drift_threshold",
      "lr_alpha",
      "louvain_iters",
      "full_snn",
      "pruning",
      "snn_similarity"
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  res <- checkKnnParams(x[names(x) %in% KNN_PARAM_NAMES])
  if (!isTRUE(res)) {
    return(res)
  }

  res <- apply_qtest_rules(
    x,
    list(
      kmeans_iters = "I1[1,)",
      batch_size = "I1[1,)",
      louvain_iters = "I1[1,)",
      drift_threshold = "N1",
      lr_alpha = "N1",
      full_snn = "B1",
      pruning = c("0", "N1")
    ),
    label = "fast clustering params",
    hint = paste(
      "kmeans_iters, batch_size and louvain_iters must be integers >= 1;",
      "drift_threshold and lr_alpha must be single numerics;",
      "full_snn must be a boolean;",
      "pruning must be NULL or a single numeric."
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  apply_choice_rules(
    x,
    list(snn_similarity = c("jaccard", "rank")),
    label = "fast clustering params"
  )
}

#' Assert SC fast clustering parameters
#'
#' @description Checkmate extension for asserting the fast clustering
#' parameters.
#'
#' @inheritParams checkScFastCluster
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#'   to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#'   [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertScFastCluster <- checkmate::makeAssertionFunction(checkScFastCluster)

#### sc type -------------------------------------------------------------------

#' Check cell marker list
#'
#' @description Checkmate extension for checking a cell marker list as returned
#' by [prepare_cell_markers()].
#'
#' @param x The list to check.
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkCellMarkerList <- function(x) {
  res <- checkmate::checkList(x, names = "named", min.len = 1)
  if (!isTRUE(res)) {
    return(res)
  }

  for (nm in names(x)) {
    entry <- x[[nm]]

    res <- checkmate::checkList(entry, names = "named", len = 3)
    if (!isTRUE(res)) {
      return(sprintf("Entry '%s': %s", nm, res))
    }

    res <- checkmate::checkNames(
      names(entry),
      must.include = c("cell_type", "positive_indices", "negative_indices")
    )
    if (!isTRUE(res)) {
      return(sprintf("Entry '%s': %s", nm, res))
    }

    res <- checkmate::checkString(entry$cell_type, min.chars = 1)
    if (!isTRUE(res)) {
      return(sprintf("Entry '%s' cell_type: %s", nm, res))
    }

    res <- checkmate::checkIntegerish(entry$positive_indices, min.len = 1)
    if (!isTRUE(res)) {
      return(sprintf("Entry '%s' positive_indices: %s", nm, res))
    }

    if (!is.null(entry$negative_indices)) {
      res <- checkmate::checkIntegerish(entry$negative_indices, min.len = 1)
      if (!isTRUE(res)) {
        return(sprintf("Entry '%s' negative_indices: %s", nm, res))
      }
    }
  }

  TRUE
}

#' Assert cell marker list
#'
#' @description Checkmate extension for asserting a cell marker list as returned
#' by [prepare_cell_markers()].
#'
#' @param x The list to assert.
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns \code{x} if the assertion is successful.
#'
#' @keywords internal
assertCellMarkerList <- checkmate::makeAssertionFunction(checkCellMarkerList)

#' Check per-cell ScType parameters
#'
#' @description Checkmate extension for checking the per-cell ScType parameters
#' as returned by [params_sctype_cells()].
#'
#' @param x The list to check.
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkSctypeCellParams <- function(x) {
  res <- checkmate::checkList(x)
  if (!isTRUE(res)) {
    return(res)
  }

  rules <- list(
    "alpha" = "N1[0,1]",
    "iterations" = "I1[0,)",
    "tolerance" = "N1(0,)",
    "calibration" = "S1",
    "score_floor" = "N1[0,)",
    "purity_threshold" = "N1[0,1]"
  )

  res <- checkmate::checkNames(names(x), must.include = names(rules))
  if (!isTRUE(res)) {
    return(res)
  }

  res <- purrr::imap_lgl(x, \(x, name) checkmate::qtest(x, rules[[name]]))
  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(sprintf("Element `%s` does not conform.", broken_elem))
  }

  res <- checkmate::checkChoice(x$calibration, c("none", "column_z"))
  if (!isTRUE(res)) {
    return(sprintf("Element `calibration`: %s", res))
  }

  TRUE
}

#' Assert per-cell ScType parameters
#'
#' @description Checkmate extension for asserting the per-cell ScType parameters
#' as returned by [params_sctype_cells()].
#'
#' @param x The list to assert.
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns \code{x} if the assertion is successful.
#'
#' @keywords internal
assertSctypeCellParams <- checkmate::makeAssertionFunction(
  checkSctypeCellParams
)

#### meld ----------------------------------------------------------------------

#' Check MELD parameters
#'
#' @description Checkmate extension for checking MELD parameters.
#'
#' @param x The list to check/assert.
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkMeldParams <- function(x) {
  res <- check_list_shape(
    x,
    c(
      "beta",
      "offset",
      "order",
      "filter",
      "chebyshev_order",
      "lap_type",
      "normalise_indicators"
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  res <- checkKnnParams(x[names(x) %in% KNN_PARAM_NAMES])
  if (!isTRUE(res)) {
    return(res)
  }

  res <- apply_qtest_rules(
    x,
    list(
      beta = "N1(0,)",
      offset = "N1[0,1]",
      order = "N1(0,)",
      chebyshev_order = "I1[2,)",
      normalise_indicators = "B1"
    ),
    label = "MELD params",
    hint = paste(
      "beta and order must be positive numerics; offset must be in [0, 1];",
      "chebyshev_order must be an integer >= 2;",
      "normalise_indicators must be a single logical."
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  apply_choice_rules(
    x,
    list(
      filter = c("heat", "laplacian"),
      lap_type = c("combinatorial", "normalised")
    ),
    label = "MELD params"
  )
}

#' Assert MELD parameters
#'
#' @description Checkmate extension for asserting MELD parameters.
#'
#' @inheritParams checkMeldParams
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertMeldParams <- checkmate::makeAssertionFunction(checkMeldParams)

### single cells (multi modal) -------------------------------------------------

#### dsb count normalisation ---------------------------------------------------

#' Check DSB parameters
#'
#' @description Checkmate extension for checking DSB parameters.
#'
#' @param x The list to check/assert.
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkScDsbParams <- function(x) {
  res <- check_list_shape(
    x,
    c(
      "denoise_counts",
      "use_isotype_controls",
      "pseudocount",
      "quantile_low",
      "quantile_high"
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  res <- apply_qtest_rules(
    x,
    list(
      denoise_counts = "B1",
      use_isotype_controls = "B1",
      pseudocount = "N1(0,)",
      quantile_low = c("N1[0,1)", "0"),
      quantile_high = c("N1(0,1]", "0")
    ),
    label = "DSB params",
    hint = paste(
      "denoise_counts and use_isotype_controls must be single logicals;",
      "pseudocount must be a single numeric > 0;",
      "quantile_low must be NULL or in [0, 1);",
      "quantile_high must be NULL or in (0, 1]."
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  if (xor(is.null(x$quantile_low), is.null(x$quantile_high))) {
    return("quantile_low and quantile_high must both be provided or both NULL.")
  }
  if (
    !is.null(x$quantile_low) &&
      !is.null(x$quantile_high) &&
      x$quantile_low >= x$quantile_high
  ) {
    return("quantile_low must be strictly less than quantile_high.")
  }

  TRUE
}

#' Assert DSB parameters
#'
#' @description Checkmate extension for asserting DSB parameters.
#'
#' @inheritParams checkScDsbParams
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertScDsbParams <- checkmate::makeAssertionFunction(checkScDsbParams)

#### wnn -----------------------------------------------------------------------

#' Check WNN parameters
#'
#' @description Checkmate extension for checking WNN parameters.
#'
#' @param x The list to check/assert.
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkScWnnParams <- function(x) {
  res <- checkmate::checkList(x)
  if (!isTRUE(res)) {
    return(res)
  }

  res <- checkKnnParams(x)
  if (!isTRUE(res)) {
    return(res)
  }

  res <- checkmate::checkNames(
    names(x),
    must.include = c(
      "k_nn",
      "knn_range",
      "sigma_method",
      "sigma_idx",
      "snn_type",
      "s_nn",
      "sd_scale",
      "kernel_power",
      "cross_const",
      "sigma_floor"
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  res <- apply_qtest_rules(
    x,
    list(
      k_nn = "I1[1,)",
      knn_range = "I1[1,)",
      sigma_idx = "I1[0,)",
      s_nn = "I1[1,)",
      sd_scale = "N1(0,)",
      kernel_power = "N1(0,)",
      cross_const = "N1[0,)",
      sigma_floor = "N1(0,)"
    ),
    label = "WNN params",
    hint = paste(
      "k_nn, knn_range and s_nn must be integers >= 1;",
      "sigma_idx must be an integer >= 0;",
      "sd_scale, kernel_power and sigma_floor must be > 0;",
      "cross_const must be >= 0."
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  apply_choice_rules(
    x,
    list(
      sigma_method = c("snn_farthest", "sigma_idx"),
      snn_type = c("full_connection", "limited")
    ),
    label = "WNN params"
  )
}

#' Assert WNN parameters
#'
#' @description Checkmate extension for asserting the WNN parameters.
#'
#' @inheritParams checkScWnnParams
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertScWnnParams <- checkmate::makeAssertionFunction(checkScWnnParams)

#### nichenet ------------------------------------------------------------------

#' Check ligand-target influence parameters
#'
#' @description Checkmate extension for checking the ligand to target
#' influence parameters.
#'
#' @param x The list to check/assert.
#'
#' @return \code{TRUE} if the check was successful, otherwise an error
#' message.
#'
#' @keywords internal
checkLigandTarget <- function(x) {
  res <- check_list_shape(
    x,
    c(
      "lr_sig_hub",
      "gr_hub",
      "ltf_cutoff",
      "damping_factor",
      "tol",
      "max_iter",
      "topology_correction",
      "secondary_targets"
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  apply_qtest_rules(
    x,
    list(
      lr_sig_hub = "N1[0,1]",
      gr_hub = "N1[0,1]",
      ltf_cutoff = "N1[0,1]",
      damping_factor = "N1[0,1]",
      tol = "N1(0,)",
      max_iter = "X1[1,)",
      topology_correction = "B1",
      secondary_targets = "B1"
    ),
    label = "ligand-target params",
    hint = paste(
      "lr_sig_hub, gr_hub, ltf_cutoff and damping_factor must be single",
      "numerics in [0, 1]; tol must be a single positive numeric;",
      "max_iter must be a single positive integer;",
      "topology_correction and secondary_targets must be single booleans."
    )
  )
}

#' Assert ligand-target influence parameters
#'
#' @inheritParams checkLigandTarget
#'
#' @param .var.name Name of the checked object to print in assertions.
#' Defaults to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is
#' successful.
#'
#' @keywords internal
assertLigandTarget <- checkmate::makeAssertionFunction(checkLigandTarget)

### spatial --------------------------------------------------------------------

#' Check that an object is a valid SpatialSample
#'
#' @description Checkmate extension for the `SpatialSample` S3 class built by
#' [bixverse::new_spatial_sample()].
#'
#' @param x The object to check.
#'
#' @return `TRUE` if the check was successful, otherwise an error message
#' string.
#'
#' @keywords internal
checkSpatialSample <- function(x) {
  res <- validate_spatial_sample(x)
  if (isTRUE(res)) TRUE else res
}

#' Assert that an object is a valid SpatialSample
#'
#' @description Checkmate assertion for the `SpatialSample` S3 class built by
#' [bixverse::new_spatial_sample()].
#'
#' @inheritParams checkSpatialSample
#'
#' @param .var.name Name of the checked object to print in assertions.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is
#' successful.
#'
#' @keywords internal
assertSpatialSample <- checkmate::makeAssertionFunction(checkSpatialSample)

#' Check Visium ingest parameters
#'
#' @description Checkmate extension for checking the Visium ingest parameters
#' built by [bixverse::params_sp_visium_io()].
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkSpVisiumIo <- function(x) {
  res <- check_list_shape(
    x,
    c("in_tissue_only", "matrix_type", "slide_file")
  )
  if (!isTRUE(res)) {
    return(res)
  }

  res <- apply_qtest_rules(
    x,
    list(
      in_tissue_only = "B1",
      matrix_type = "S1",
      slide_file = c("0", "S1")
    ),
    label = "Visium IO params",
    hint = paste(
      "in_tissue_only must be a single boolean, matrix_type a single string",
      "and slide_file either NULL or a single string."
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  apply_choice_rules(
    x,
    list(matrix_type = c("auto", "raw", "filtered")),
    label = "Visium IO params",
    hint = "matrix_type must be one of 'auto', 'raw' or 'filtered'."
  )
}

#' Assert Visium ingest parameters
#'
#' @description Checkmate assertion for the Visium ingest parameters built by
#' [bixverse::params_sp_visium_io()].
#'
#' @inheritParams checkSpVisiumIo
#'
#' @param .var.name Name of the checked object to print in assertions.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is
#' successful.
#'
#' @keywords internal
assertSpVisiumIo <- checkmate::makeAssertionFunction(checkSpVisiumIo)

#' Check spatial graph parameters
#'
#' @description Checkmate extension for checking the spatial graph parameters
#' built by [bixverse::params_sp_graph()].
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkSpGraph <- function(x) {
  res <- check_list_shape(
    x,
    c(
      "layout",
      "weighting",
      "connectivity",
      "k",
      "radius",
      "power",
      "bandwidth",
      "row_normalise"
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  res <- apply_qtest_rules(
    x,
    list(
      layout = "S1",
      weighting = "S1",
      connectivity = "I1",
      k = "I1[1,)",
      radius = c("0", "N1(0,)"),
      power = "N1(0,)",
      bandwidth = c("0", "N1(0,)"),
      row_normalise = "B1"
    ),
    label = "spatial graph params",
    hint = paste(
      "layout and weighting must be single strings, connectivity and k single",
      "integers, power a positive float, radius and bandwidth either NULL or a",
      "positive float and row_normalise a single boolean."
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  res <- apply_choice_rules(
    x,
    list(
      layout = c("hex", "square", "knn", "radius"),
      weighting = c("binary", "inverse_distance", "gaussian"),
      connectivity = c(4L, 8L)
    ),
    label = "spatial graph params",
    hint = paste(
      "layout must be one of 'hex', 'square', 'knn' or 'radius', weighting one",
      "of 'binary', 'inverse_distance' or 'gaussian' and connectivity 4 or 8."
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  if (identical(x$layout, "radius") && is.null(x$radius)) {
    return("The layout 'radius' in spatial graph params needs a `radius`.")
  }

  return(TRUE)
}

#' Assert spatial graph parameters
#'
#' @description Checkmate assertion for the spatial graph parameters built by
#' [bixverse::params_sp_graph()].
#'
#' @inheritParams checkSpGraph
#'
#' @param .var.name Name of the checked object to print in assertions.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is
#' successful.
#'
#' @keywords internal
assertSpGraph <- checkmate::makeAssertionFunction(checkSpGraph)

#' Check Moran's I parameters
#'
#' @description Checkmate extension for checking the Moran's I parameters built
#' by [bixverse::params_sp_svg()].
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkSpSvg <- function(x) {
  res <- check_list_shape(x, "assay")
  if (!isTRUE(res)) {
    return(res)
  }

  res <- apply_qtest_rules(
    x,
    list(assay = "S1"),
    label = "Moran's I params",
    hint = "assay must be a single string."
  )
  if (!isTRUE(res)) {
    return(res)
  }

  apply_choice_rules(
    x,
    list(assay = c("raw", "norm")),
    label = "Moran's I params",
    hint = "assay must be either 'raw' or 'norm'."
  )
}

#' Assert Moran's I parameters
#'
#' @description Checkmate assertion for the Moran's I parameters built by
#' [bixverse::params_sp_svg()].
#'
#' @inheritParams checkSpSvg
#'
#' @param .var.name Name of the checked object to print in assertions.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is
#' successful.
#'
#' @keywords internal
assertSpSvg <- checkmate::makeAssertionFunction(checkSpSvg)

#' Check SPARK-X parameters
#'
#' @description Checkmate extension for checking the SPARK-X parameters built by
#' [bixverse::params_sp_sparkx()].
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkSpSparkx <- function(x) {
  res <- check_list_shape(
    x,
    c("kernels", "n_landmarks", "bandwidth_subsample")
  )
  if (!isTRUE(res)) {
    return(res)
  }

  res <- apply_qtest_rules(
    x,
    list(
      kernels = c("0", "l+"),
      n_landmarks = "I1[1,)",
      bandwidth_subsample = "I1[1,)"
    ),
    label = "SPARK-X params",
    hint = paste(
      "kernels must be NULL or a non-empty list, n_landmarks and",
      "bandwidth_subsample positive integers."
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  if (is.null(x$kernels)) {
    return(TRUE)
  }

  ok <- purrr::map_lgl(x$kernels, \(kernel) {
    checkmate::testList(kernel) &&
      all(c("kernel", "bandwidth") %in% names(kernel)) &&
      checkmate::testChoice(kernel$kernel, c("gaussian", "cosine")) &&
      checkmate::qtest(kernel$bandwidth, "N1(0,)")
  })
  if (!all(ok)) {
    return(sprintf(
      paste(
        "Kernel %i in SPARK-X params is invalid. Every kernel needs a `kernel`",
        "of 'gaussian' or 'cosine' and a positive `bandwidth`."
      ),
      which(!ok)[1]
    ))
  }

  return(TRUE)
}

#' Assert SPARK-X parameters
#'
#' @description Checkmate assertion for the SPARK-X parameters built by
#' [bixverse::params_sp_sparkx()].
#'
#' @inheritParams checkSpSparkx
#'
#' @param .var.name Name of the checked object to print in assertions.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is
#' successful.
#'
#' @keywords internal
assertSpSparkx <- checkmate::makeAssertionFunction(checkSpSparkx)

#' Check neighbourhood enrichment parameters
#'
#' @description Checkmate extension for checking the neighbourhood enrichment
#' parameters built by [bixverse::params_sp_nhood()].
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkSpNhood <- function(x) {
  res <- check_list_shape(x, c("n_perm", "symmetrise"))
  if (!isTRUE(res)) {
    return(res)
  }

  apply_qtest_rules(
    x,
    list(n_perm = "I1[1,)", symmetrise = "B1"),
    label = "neighbourhood enrichment params",
    hint = "n_perm must be a positive integer and symmetrise a single boolean."
  )
}

#' Assert neighbourhood enrichment parameters
#'
#' @description Checkmate assertion for the neighbourhood enrichment parameters
#' built by [bixverse::params_sp_nhood()].
#'
#' @inheritParams checkSpNhood
#'
#' @param .var.name Name of the checked object to print in assertions.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is
#' successful.
#'
#' @keywords internal
assertSpNhood <- checkmate::makeAssertionFunction(checkSpNhood)

#' Check histology image feature parameters
#'
#' @description Checkmate extension for checking the image feature parameters
#' built by [bixverse::params_sp_image()].
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkSpImage <- function(x) {
  res <- check_list_shape(
    x,
    c(
      "tile_scale",
      "glcm_levels",
      "glcm_offsets_dy",
      "glcm_offsets_dx",
      "stain_haem",
      "stain_eosin"
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  res <- apply_qtest_rules(
    x,
    list(
      tile_scale = "N1(0,)",
      glcm_levels = "I1[2,255]",
      glcm_offsets_dy = c("0", "I+"),
      glcm_offsets_dx = c("0", "I+"),
      stain_haem = c("0", "N3"),
      stain_eosin = c("0", "N3")
    ),
    label = "image feature params",
    hint = paste(
      "tile_scale must be a positive float, glcm_levels an integer between 2",
      "and 255, the GLCM offsets NULL or integer vectors and the stain vectors",
      "NULL or numerics of length three."
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  if (is.null(x$glcm_offsets_dy) != is.null(x$glcm_offsets_dx)) {
    return(paste(
      "`glcm_offsets_dy` and `glcm_offsets_dx` in image feature params must be",
      "given together."
    ))
  }
  if (
    !is.null(x$glcm_offsets_dy) &&
      length(x$glcm_offsets_dy) != length(x$glcm_offsets_dx)
  ) {
    return(paste(
      "`glcm_offsets_dy` and `glcm_offsets_dx` in image feature params must",
      "have the same length."
    ))
  }
  if (is.null(x$stain_haem) != is.null(x$stain_eosin)) {
    return(paste(
      "`stain_haem` and `stain_eosin` in image feature params must be given",
      "together."
    ))
  }

  return(TRUE)
}

#' Assert histology image feature parameters
#'
#' @description Checkmate assertion for the image feature parameters built by
#' [bixverse::params_sp_image()].
#'
#' @inheritParams checkSpImage
#'
#' @param .var.name Name of the checked object to print in assertions.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is
#' successful.
#'
#' @keywords internal
assertSpImage <- checkmate::makeAssertionFunction(checkSpImage)
