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
  res <- purrr::keep(
    res,
    ~ {
      !is.logical(.x)
    }
  )
  if (length(res) == 0) {
    return(TRUE)
  } else {
    return(res[[1]])
  }
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
  # Checkmate extension
  res <- checkmate::checkList(x)
  if (!isTRUE(res)) {
    return(res)
  }
  res <- checkmate::checkNames(
    names(x),
    must.include = c("epsilon", "min_cor", "fdr_threshold", "verbose")
  )
  if (!isTRUE(res)) {
    return(res)
  }
  rules <- list(
    "epsilon" = "R1",
    "min_cor" = "R1[0, 1]",
    "fdr_threshold" = "R1[0, 1]",
    "verbose" = "B1"
  )
  res <- purrr::imap_lgl(x, \(x, name) {
    checkmate::qtest(x, rules[[name]])
  })
  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(
      sprintf(
        paste(
          "The following element `%s` in graph params does not conform to the",
          "expected format. min_cor and fdr_threshold need to be between 0 and",
          "1, epsilon a double and .verbose a boolean."
        ),
        broken_elem
      )
    )
  }
  return(TRUE)
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
  # Checkmate extension
  res <- checkmate::checkList(x)
  if (!isTRUE(res)) {
    return(res)
  }
  res <- checkmate::checkNames(
    names(x),
    must.include = c("min_res", "max_res", "number_res")
  )
  if (!isTRUE(res)) {
    return(res)
  }
  rules <- list("min_res" = "R1", "max_res" = "R1", "number_res" = "I1")
  res <- purrr::imap_lgl(x, \(x, name) {
    checkmate::qtest(x, rules[[name]])
  })
  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(
      sprintf(
        paste(
          "The following element `%s` in resolution params does not conform to",
          "the expected format. min_res and max_res need to be doubles and",
          "number res needs to be an integer."
        ),
        broken_elem
      )
    )
  }
  return(TRUE)
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
  # Checkmate extension
  res <- checkmate::checkList(x)
  if (!isTRUE(res)) {
    return(res)
  }
  res <- checkmate::checkNames(
    names(x),
    must.include = c("maxit", "alpha", "max_tol", "verbose")
  )
  if (!isTRUE(res)) {
    return(res)
  }
  rules <- list(
    "maxit" = "I1",
    "alpha" = "R1[1, 2]",
    "max_tol" = "R1(0, 1)",
    "verbose" = "B1"
  )
  res <- purrr::imap_lgl(x, \(x, name) {
    checkmate::qtest(x, rules[[name]])
  })
  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(
      sprintf(
        paste(
          "The following element `%s` in resolution params does not conform to",
          "the expected format. maxit needs to be an integer, alpha between 1",
          "and 2, 0 < max_tol < 1, and verbose a boolean."
        ),
        broken_elem
      )
    )
  }
  return(TRUE)
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
  # Checkmate extension
  res <- checkmate::checkList(x)
  if (!isTRUE(res)) {
    return(res)
  }
  res <- checkmate::checkNames(
    names(x),
    must.include = c("max_no_comp", "steps", "custom_seq")
  )
  if (!isTRUE(res)) {
    return(res)
  }
  rules <- list(
    "max_no_comp" = "I1",
    "steps" = "I1",
    "custom_seq" = c("0", "I+")
  )
  res <- purrr::imap_lgl(x, \(x, name) {
    checkmate::qtest(x, rules[[name]])
  })
  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(
      sprintf(
        paste(
          "The following element `%s` in resolution params does not conform to",
          "the expected format. max_no_comp and steps need to be integers, and",
          "custom sequence either NULL or a vector of integers."
        ),
        broken_elem
      )
    )
  }
  return(TRUE)
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
  # Checkmate extension
  res <- checkmate::checkList(x)
  if (!isTRUE(res)) {
    return(res)
  }
  res <- checkmate::checkNames(
    names(x),
    must.include = c("cross_validate", "random_init", "folds")
  )
  if (!isTRUE(res)) {
    return(res)
  }
  rules <- list(
    "cross_validate" = "B1",
    "random_init" = "I1",
    "folds" = "I1"
  )
  res <- purrr::imap_lgl(x, \(x, name) {
    checkmate::qtest(x, rules[[name]])
  })
  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(
      sprintf(
        paste(
          "The following element `%s` in resolution params does not conform to",
          "the expected format. random_init and steps folds need to be",
          "integers, and cross_validate a boolean."
        ),
        broken_elem
      )
    )
  }
  return(TRUE)
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
  # Checkmate extension
  res <- checkmate::checkList(x)
  if (!isTRUE(res)) {
    return(res)
  }
  res <- checkmate::checkNames(
    names(x),
    must.include = c(
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
  res <- checkmate::checkChoice(
    x[['threshold_type']],
    c("prop_based", "pval_based")
  )
  if (!isTRUE(res)) {
    return(res)
  }
  q_rules <- list(
    "max_nodes" = sprintf("I1[%i,)", x$min_nodes),
    "min_nodes" = "I1",
    "min_seed_nodes" = "I1",
    "initial_res" = "N1",
    "network_threshold" = "N1(0, 1]",
    "pval_threshold" = "N1(0, 1]"
  )
  res <- purrr::imap_lgl(x, \(x, name) {
    if (name %in% names(q_rules)) {
      checkmate::qtest(x, q_rules[[name]])
    } else {
      TRUE
    }
  })

  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(
      sprintf(
        paste(
          "The following element `%s` in community params does not conform to",
          "the expected format. min_nodes, max_nodes and min_seed_genes need to",
          "be integers (with max_nodes > min_nodes), initial resolution a",
          "double, and network_threshold and pval_threshold doubles between 0 and 1."
        ),
        broken_elem
      )
    )
  }
  return(TRUE)
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
  # Checkmate extension
  res <- checkmate::checkList(x)
  if (!isTRUE(res)) {
    return(res)
  }
  res <- checkmate::checkNames(
    names(x),
    must.include = c("min_size", "max_size", "gsea_param", "sample_size", "eps")
  )
  if (!isTRUE(res)) {
    return(res)
  }
  rules <- list(
    "min_size" = "I1[3,)",
    "max_size" = "I1[4,)",
    "gsea_param" = "N1",
    "sample_size" = "I1",
    "eps" = "N1"
  )
  res <- purrr::imap_lgl(x, \(x, name) {
    checkmate::qtest(x, rules[[name]])
  })
  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(
      sprintf(
        paste(
          "The following element `%s` in GSEA params does not conform to the",
          "expected format. min_size and max_size need to be integers (with",
          "max_size > min_size and min_size >= 3L),",
          "gsea_param being a double, sample_size an integer and eps a float."
        ),
        broken_elem
      )
    )
  }
  return(TRUE)
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
  # Checkmate extension
  res <- checkmate::checkList(x)
  if (!isTRUE(res)) {
    return(res)
  }
  res <- checkmate::checkNames(
    names(x),
    must.include = c(
      "tau",
      "min_size",
      "max_size",
      "max_diff",
      "abs_rank"
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }
  rules <- list(
    "tau" = "N1",
    "min_size" = "I1[3,)",
    "max_size" = "I1[4,)",
    "max_diff" = "B1",
    "abs_rank" = "B1"
  )
  res <- purrr::imap_lgl(x, \(x, name) {
    checkmate::qtest(x, rules[[name]])
  })
  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(
      sprintf(
        paste(
          "The following element `%s` in GSVA params does not conform to the",
          "expected format. min_size and max_size need to be integers (with",
          "max_size > min_size and min_size >= 3L),",
          "tau being a double, max_diff and abs_rank booleans."
        ),
        broken_elem
      )
    )
  }
  return(TRUE)
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
  # Checkmate extension
  res <- checkmate::checkList(x)
  if (!isTRUE(res)) {
    return(res)
  }
  res <- checkmate::checkNames(
    names(x),
    must.include = c(
      "alpha",
      "min_size",
      "max_size",
      "normalise"
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }
  rules <- list(
    "alpha" = "N1(0,1)",
    "min_size" = "I1[3,)",
    "max_size" = "I1[4,)",
    "normalise" = "B1"
  )
  res <- purrr::imap_lgl(x, \(x, name) {
    checkmate::qtest(x, rules[[name]])
  })
  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(
      sprintf(
        paste(
          "The following element `%s` in ssGSEA params does not conform to the",
          "expected format. min_size and max_size need to be integers (with",
          "max_size > min_size and min_size >= 3L),",
          "alpha being a double (between 0 and 1), and normalise a boolean."
        ),
        broken_elem
      )
    )
  }
  return(TRUE)
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
  # Checkmate extension
  res <- checkmate::checkList(x)
  if (!isTRUE(res)) {
    return(res)
  }
  res <- checkmate::checkNames(
    names(x),
    must.include = c(
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
  # qtest checks
  qtest_rules <- list(
    epsilon = "N1",
    k_min = "I1",
    k_max = "I1",
    junk_module_threshold = "N1",
    min_size = c("I1", "0")
  )
  q_test_res <- purrr::imap_lgl(x, \(x, name) {
    if (name %in% names(qtest_rules)) {
      checkmate::qtest(x, qtest_rules[[name]])
    } else {
      TRUE
    }
  })
  if (!isTRUE(all(q_test_res))) {
    broken_elem <- names(q_test_res)[which(!q_test_res)][1]
    return(
      sprintf(
        paste(
          "The following element `%s` in CoReMo params does not conform to the",
          "expected format k_min and k_max need to be integers, min_size an",
          "integer or NULL, junk_module_threshold a float and epsilon a float."
        ),
        broken_elem
      )
    )
  }
  # test
  test_choice_rules <- list(
    rbf_func = c("gaussian", "inverse_quadratic", "bump"),
    cor_method = c("spearman", "pearson")
  )
  test_choice_res <- purrr::imap_lgl(x, \(x, name) {
    if (name %in% names(test_choice_rules)) {
      checkmate::testChoice(x, test_choice_rules[[name]])
    } else {
      TRUE
    }
  })
  if (!isTRUE(all(test_choice_res))) {
    broken_elem <- names(test_choice_res)[which(!test_choice_res)][1]
    return(
      sprintf(
        paste0(
          "The following element `%s` in CoReMo params does not use one of the",
          "expected choices. Please double check the documentation."
        ),
        broken_elem
      )
    )
  }
  return(TRUE)
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
  # Checkmate extension
  res <- checkmate::checkList(x)
  if (!isTRUE(res)) {
    return(res)
  }
  res <- checkmate::checkNames(
    names(x),
    must.include = c(
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
  rules <- list(
    "sparsity" = "I1",
    "dict_size" = "I1",
    "alpha" = "N1",
    "beta" = "N1",
    "max_iter" = "I1",
    "k_neighbours" = "I1",
    "admm_iter" = "I1",
    "rho" = "N1"
  )
  res <- purrr::imap_lgl(x, \(x, name) {
    checkmate::qtest(x, rules[[name]])
  })
  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(
      sprintf(
        paste(
          "The following element `%s` in DGRDL params does not conform to the",
          "expected format. sparsity, dict_size, max_iter, k_neighbours, and",
          "admm_iter are expected to be integers; alpha, beta, rho are",
          "expected to be floats."
        ),
        broken_elem
      )
    )
  }
  return(TRUE)
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
  # Check it's a list
  res <- checkmate::checkList(x)
  if (!isTRUE(res)) {
    return(res)
  }

  # Check required names
  res <- checkmate::checkNames(
    names(x),
    must.include = c(
      "k",
      "t",
      "mu",
      "alpha",
      "normalise",
      "distance_metric"
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  # qtest checks
  qtest_rules <- list(
    k = "I1",
    t = "I1",
    mu = "N1[0,1]",
    alpha = "N1",
    normalise = "B1"
  )

  q_test_res <- purrr::imap_lgl(x, \(x, name) {
    if (name %in% names(qtest_rules)) {
      checkmate::qtest(x, qtest_rules[[name]])
    } else {
      TRUE
    }
  })

  if (!isTRUE(all(q_test_res))) {
    broken_elem <- names(q_test_res)[which(!q_test_res)][1]
    return(
      sprintf(
        paste(
          "The following element `%s` in SNF params does not conform to the",
          "expected format. k and t need to be positive integers, mu a float",
          "in [0, 1], alpha a float, and normalise a boolean."
        ),
        broken_elem
      )
    )
  }

  # test choice rules
  test_choice_rules <- list(
    distance_metric = c("euclidean", "manhattan", "canberra", "cosine")
  )

  test_choice_res <- purrr::imap_lgl(x, \(x, name) {
    if (name %in% names(test_choice_rules)) {
      checkmate::testChoice(x, test_choice_rules[[name]])
    } else {
      TRUE
    }
  })

  if (!isTRUE(all(test_choice_res))) {
    broken_elem <- names(test_choice_res)[which(!test_choice_res)][1]
    return(
      sprintf(
        paste0(
          "The following element `%s` in SNF params does not use one of the",
          " expected choices. Please double check the documentation."
        ),
        broken_elem
      )
    )
  }

  return(TRUE)
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
  # Check it's a list
  res <- checkmate::checkList(x)
  if (!isTRUE(res)) {
    return(res)
  }

  # Check required names
  res <- checkmate::checkNames(
    names(x),
    must.include = c(
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

  # Validate types
  rules <- list(
    "auc_threshold" = "N1[0,1]",
    "nes_threshold" = "N1",
    "rcc_method" = "S1",
    "high_conf_cats" = "S+",
    "low_conf_cats" = "S+"
  )

  res <- purrr::imap_lgl(x, \(val, name) {
    checkmate::qtest(val, rules[[name]])
  })

  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(
      sprintf(
        paste(
          "The following element `%s` in CisTarget params does not conform to",
          "the expected format. auc_threshold must be numeric [0,1];",
          "nes_threshold must be numeric; rcc_method must be a single string;",
          "high_conf_cats and low_conf_cats must be character vectors."
        ),
        broken_elem
      )
    )
  }

  if (!checkmate::testChoice(x$rcc_method, c("approx", "icistarget"))) {
    return("rcc_method must be either 'approx' or 'icistarget'")
  }

  return(TRUE)
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
  res <- checkmate::checkList(x)
  if (!isTRUE(res)) {
    return(res)
  }

  res <- checkmate::checkNames(
    names(x),
    must.include = c(
      "alpha",
      "iter",
      "tolerance",
      "symmetrise",
      "symmetry_strategy"
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  rules <- list(
    "alpha" = "R1[0, 1]",
    "iter" = "I1[1,]",
    "tolerance" = "R1",
    "symmetrise" = "B1",
    "symmetry_strategy" = "S1"
  )
  res <- purrr::imap_lgl(x[names(rules)], \(x, name) {
    checkmate::qtest(x, rules[[name]])
  })
  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(sprintf(
      paste(
        "The following element `%s` in label propagation params does not",
        "conform to the expected format. alpha must be in [0, 1], iter a",
        "positive integer, tolerance a double, symmetrise a boolean, and",
        "symmetry_strategy a string."
      ),
      broken_elem
    ))
  }

  if (
    !is.null(x[["max_hops"]]) && !checkmate::qtest(x[["max_hops"]], "I1[0,]")
  ) {
    return("`max_hops` must be a positive integer or NULL.")
  }

  return(TRUE)
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
  # If required_params not specified, check all that are present
  if (!is.null(required_params)) {
    res <- checkmate::checkNames(names(x), must.include = required_params)
    if (!isTRUE(res)) {
      return(res)
    }
  }

  # integer rules
  integer_rules <- list(
    "k" = "I1[0,)",
    "n_trees" = "I1[1,)",
    "search_budget" = c("0", "I1[1,)"),
    "m" = "I1[1,)",
    "ef_construction" = "I1[1,)",
    "ef_search" = "I1[1,)",
    "ef_budget" = c("0", "I1[1,)"),
    "n_list" = c("0", "I1[1,)"),
    "n_probe" = c("0", "I1[1,)")
  )

  res <- purrr::imap_lgl(x, \(val, name) {
    if (name %in% names(integer_rules)) {
      checkmate::qtest(val, integer_rules[[name]])
    } else {
      TRUE
    }
  })

  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(
      sprintf(
        paste(
          "The following kNN parameter `%s` is incorrect:",
          "k must be >= 0; n_trees, m, ef_construction, ef_search must be >= 1;",
          "search_budget, ef_budget, n_list and n_probe must be NULL or >= 1;"
        ),
        broken_elem
      )
    )
  }

  # numeric rules
  numeric_rules <- list(
    "delta" = "N1[0,1]",
    "diversify_prob" = "N1[0,1]"
  )

  res <- purrr::imap_lgl(x, \(val, name) {
    if (name %in% names(numeric_rules)) {
      checkmate::qtest(val, numeric_rules[[name]])
    } else {
      TRUE
    }
  })

  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(
      sprintf(
        "The kNN parameter `%s` must be in [0,1].",
        broken_elem
      )
    )
  }

  # choice rules
  test_choice_rules <- list(
    knn_method = c("annoy", "hnsw", "nndescent", "exhaustive", "ivf"),
    ann_dist = c("euclidean", "cosine")
  )

  test_choice_res <- purrr::imap_lgl(x, \(val, name) {
    if (name %in% names(test_choice_rules)) {
      checkmate::testChoice(val, test_choice_rules[[name]])
    } else {
      TRUE
    }
  })

  if (!isTRUE(all(test_choice_res))) {
    broken_elem <- names(test_choice_res)[which(!test_choice_res)][1]
    return(
      sprintf(
        paste(
          "The kNN parameter `%s` is not one of the expected choices.",
          "Please check the documentation."
        ),
        broken_elem
      )
    )
  }

  return(TRUE)
}

#### synthetic data ------------------------------------------------------------

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
  res <- checkmate::checkList(x)
  if (!isTRUE(res)) {
    return(res)
  }

  res <- checkmate::checkNames(
    names(x),
    must.include = c(
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

  rules <- list(
    "n_cells" = "I1",
    "n_genes" = "I1",
    "n_batches" = "I1",
    "n_samples" = c("0", "I1")
  )

  res <- purrr::imap_lgl(x, \(x, name) {
    if (name %in% names(rules)) {
      checkmate::qtest(x, rules[[name]])
    } else {
      TRUE
    }
  })

  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(
      sprintf(
        paste(
          "The following element `%s` in synthetic data is incorrect:",
          "n_cells, n_genes and n_batches need to be integers. n_samples an",
          "integer or NULL."
        ),
        broken_elem
      )
    )
  }

  res <- checkmate::checkList(
    x$marker_genes,
    types = "list",
    names = "named"
  )
  if (!isTRUE(res)) {
    return("marker_genes must be a named list of lists.")
  }

  res <- checkmate::checkChoice(
    x[["batch_effect_strength"]],
    c("strong", "medium", "weak")
  )
  if (!isTRUE(res)) {
    return(res)
  }

  res <- if (!is.null(x[["sample_bias"]])) {
    checkmate::checkChoice(
      x[["sample_bias"]],
      c("even", "slightly_uneven", "very_uneven")
    )
  } else {
    TRUE
  }
  if (!isTRUE(res)) {
    return(res)
  }

  return(TRUE)
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
  res <- checkmate::checkList(x)
  if (!isTRUE(res)) {
    return(res)
  }
  res <- checkmate::checkNames(
    names(x),
    must.include = c(
      "path_mtx",
      "path_obs",
      "path_var",
      "cells_as_rows",
      "has_hdr"
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }
  res <- purrr::map_lgl(c("path_mtx", "path_obs", "path_var"), \(n) {
    checkmate::testFileExists(x[[n]])
  })
  if (!isTRUE(all(res))) {
    return(paste(
      "Some of the files specified in the config for mtx ingest are not",
      "existing. Please check the provided params."
    ))
  }
  res <- checkmate::qtest(x[["cells_as_rows"]], "B1")
  if (!isTRUE(res)) {
    return(res)
  }
  res <- checkmate::qtest(x[["has_hdr"]], "B1")
  if (!isTRUE(res)) {
    return(res)
  }
  return(TRUE)
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
  # Checkmate extension
  res <- checkmate::checkList(x)
  if (!isTRUE(res)) {
    return(res)
  }
  res <- checkmate::checkNames(
    names(x),
    must.include = c(
      "min_unique_genes",
      "min_lib_size",
      "min_cells",
      "target_size"
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }
  rules <- list(
    "min_unique_genes" = "I1",
    "min_lib_size" = "I1",
    "min_cells" = "I1",
    "target_size" = "N1"
  )
  res <- purrr::imap_lgl(x, \(x, name) {
    checkmate::qtest(x, rules[[name]])
  })
  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(
      sprintf(
        paste(
          "The following element `%s` in min single cell QC params does not 
          conform to the expected format. min_unique_genes, min_lib_size, ",
          "min_cells need to be integers and target_size needs to be float."
        ),
        broken_elem
      )
    )
  }
  return(TRUE)
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
  res <- checkmate::checkList(x)
  if (!isTRUE(res)) {
    return(res)
  }

  res <- checkmate::checkNames(
    names(x),
    must.include = c(
      "log_transform",
      "mean_center",
      "normalise_variance",
      "target_size",
      "min_gene_var_pctl",
      "hvg_method",
      "loess_span",
      "clip_max",
      "sim_doublet_ratio",
      "expected_doublet_rate",
      "stdev_doublet_rate",
      "manual_threshold",
      "n_bins",
      "no_pcs",
      "random_svd"
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  # Check kNN parameters
  knn_params <- x[names(x) %in% KNN_PARAM_NAMES]
  res <- checkKnnParams(knn_params)
  if (!isTRUE(res)) {
    return(res)
  }

  # Integer rules (non-kNN)
  integer_rules <- list(
    "no_pcs" = "I1[1,)",
    "n_bins" = "I1[10,)"
  )

  res <- purrr::imap_lgl(x, \(x, name) {
    if (name %in% names(integer_rules)) {
      checkmate::qtest(x, integer_rules[[name]])
    } else {
      TRUE
    }
  })

  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(
      sprintf(
        paste(
          "The following element `%s` in Scrublet parameters is incorrect:",
          "no_pcs must be >= 1; n_bins must be >= 10."
        ),
        broken_elem
      )
    )
  }

  # Numeric rules (non-kNN)
  numeric_rules <- list(
    "min_gene_var_pctl" = "N1[0,1]",
    "loess_span" = "N1(0,)",
    "sim_doublet_ratio" = "N1(0,)",
    "expected_doublet_rate" = "N1[0,1]",
    "stdev_doublet_rate" = "N1[0,1]",
    "target_size" = "N1[0,)"
  )

  res <- purrr::imap_lgl(x, \(x, name) {
    if (name %in% names(numeric_rules)) {
      checkmate::qtest(x, numeric_rules[[name]])
    } else {
      TRUE
    }
  })

  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(
      sprintf(
        paste(
          "The following element `%s` in Scrublet parameters is incorrect:",
          "min_gene_var_pctl, expected_doublet_rate and stdev_doublet_rate",
          "must be in [0,1]; loess_span and sim_doublet_ratio must be > 0;",
          "target_size must be a numeric >= 1"
        ),
        broken_elem
      )
    )
  }

  # Boolean rules
  boolean_rules <- c(
    "random_svd",
    "log_transform",
    "mean_center",
    "normalise_variance"
  )

  res <- purrr::map_lgl(boolean_rules, \(name) {
    checkmate::qtest(x[[name]], "B1")
  })

  if (!isTRUE(all(res))) {
    broken_elem <- boolean_rules[which(!res)][1]
    return(
      sprintf(
        "The element `%s` in Scrublet parameters must be a boolean (TRUE/FALSE).",
        broken_elem
      )
    )
  }

  # Optional numeric rules (can be NULL)
  optional_rules <- list(
    "clip_max" = c("0", "N1(0,)"),
    "manual_threshold" = c("0", "N1[0,)")
  )

  res <- purrr::imap_lgl(x, \(x, name) {
    if (name %in% names(optional_rules)) {
      checkmate::qtest(x, optional_rules[[name]])
    } else {
      TRUE
    }
  })

  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(
      sprintf(
        paste(
          "The following element `%s` in Scrublet parameters is incorrect:",
          "clip_max and manual_threshold must be NULL or positive numeric",
          "values."
        ),
        broken_elem
      )
    )
  }

  # Choice rules
  res <- checkmate::testChoice(x[["hvg_method"]], c("vst", "mvb", "dispersion"))
  if (!isTRUE(res)) {
    return("hvg_method must be one of: vst, mvb, dispersion.")
  }

  return(TRUE)
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
  res <- checkmate::checkList(x)
  if (!isTRUE(res)) {
    return(res)
  }

  res <- checkmate::checkNames(
    names(x),
    must.include = c(
      "log_transform",
      "mean_center",
      "normalise_variance",
      "target_size",
      "min_gene_var_pctl",
      "hvg_method",
      "loess_span",
      "clip_max",
      "boost_rate",
      "replace",
      "no_pcs",
      "random_svd",
      "resolution",
      "n_iters",
      "p_thresh",
      "voter_thresh"
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  # kNN
  knn_params <- x[names(x) %in% KNN_PARAM_NAMES]
  res <- checkKnnParams(knn_params)
  if (!isTRUE(res)) {
    return(res)
  }

  # Integer rules
  integer_rules <- list(
    "no_pcs" = "I1[1,)",
    "n_iters" = "I1[1,)"
  )

  res <- purrr::imap_lgl(x, \(x, name) {
    if (name %in% names(integer_rules)) {
      checkmate::qtest(x, integer_rules[[name]])
    } else {
      TRUE
    }
  })

  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(
      sprintf(
        paste(
          "The following element `%s` in Boost parameters is incorrect:",
          "no_pcs and n_iters must be >= 1."
        ),
        broken_elem
      )
    )
  }

  # Numeric rules
  numeric_rules <- list(
    "min_gene_var_pctl" = "N1[0,1]",
    "loess_span" = "N1(0,)",
    "boost_rate" = "N1[0,1]",
    "resolution" = "N1(0,)",
    "p_thresh" = "N1(0,)",
    "voter_thresh" = "N1[0,1]",
    "target_size" = "N1(0,)"
  )

  res <- purrr::imap_lgl(x, \(x, name) {
    if (name %in% names(numeric_rules)) {
      checkmate::qtest(x, numeric_rules[[name]])
    } else {
      TRUE
    }
  })

  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(
      sprintf(
        paste(
          "The following element `%s` in Boost parameters is incorrect:",
          "min_gene_var_pctl, boost_rate and voter_thresh must be in [0,1];",
          "loess_span, resolution and p_thresh must be > 0;",
          "target_size must be > 1."
        ),
        broken_elem
      )
    )
  }

  # Boolean rules
  boolean_rules <- c(
    "log_transform",
    "mean_center",
    "normalise_variance",
    "replace",
    "random_svd"
  )

  res <- purrr::map_lgl(boolean_rules, \(name) {
    checkmate::qtest(x[[name]], "B1")
  })

  if (!isTRUE(all(res))) {
    broken_elem <- boolean_rules[which(!res)][1]
    return(
      sprintf(
        "The element `%s` in Boost parameters must be a boolean (TRUE/FALSE).",
        broken_elem
      )
    )
  }

  # Optional numeric rules (can be NULL)
  optional_rules <- list(
    "clip_max" = c("0", "N1(0,)")
  )

  res <- purrr::imap_lgl(x, \(x, name) {
    if (name %in% names(optional_rules)) {
      checkmate::qtest(x, optional_rules[[name]])
    } else {
      TRUE
    }
  })

  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(
      sprintf(
        paste(
          "The following element `%s` in Boost parameters is incorrect:",
          "clip_max must be NULL or a positive numeric value."
        ),
        broken_elem
      )
    )
  }

  # Choice rules
  res <- checkmate::testChoice(x[["hvg_method"]], c("vst", "mvb", "dispersion"))
  if (!isTRUE(res)) {
    return("hvg_method must be one of: vst, mvb, dispersion.")
  }

  return(TRUE)
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
  res <- checkmate::checkList(x)
  if (!isTRUE(res)) {
    return(res)
  }
  res <- checkmate::checkNames(
    names(x),
    must.include = c(
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
      "n_bins"
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }
  knn_params <- x[names(x) %in% KNN_PARAM_NAMES]
  res <- checkKnnParams(knn_params)
  if (!isTRUE(res)) {
    return(res)
  }
  integer_rules <- list(
    "n_genes" = "I1[1,)",
    "no_pcs" = "I1[1,)",
    "cluster_iters" = "I1[1,)",
    "n_iterations" = "I1[1,)",
    "n_trees" = "I1[1,)",
    "max_depth" = "I1[1,)",
    "min_samples_leaf" = "I1[1,)",
    "cv_folds" = "I1[2,)",
    "cv_early_stop" = "I1[1,)",
    "include_pcs" = "I1[1,)",
    "n_bins" = "I1[10,)"
  )
  res <- purrr::imap_lgl(x, \(x, name) {
    if (name %in% names(integer_rules)) {
      checkmate::qtest(x, integer_rules[[name]])
    } else {
      TRUE
    }
  })
  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(sprintf(
      "The element `%s` in scDblFinder parameters is incorrect.",
      broken_elem
    ))
  }
  numeric_rules <- list(
    "doublet_ratio" = "N1(0,)",
    "heterotypic_bias" = "N1[0,1]",
    "cluster_resolution" = "N1(0,)",
    "learning_rate" = "N1(0,)",
    "subsample_rate" = "N1(0,1]",
    "se_fraction" = "N1[0,)",
    "target_size" = "N1(0,)"
  )
  res <- purrr::imap_lgl(x, \(x, name) {
    if (name %in% names(numeric_rules)) {
      checkmate::qtest(x, numeric_rules[[name]])
    } else {
      TRUE
    }
  })
  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(sprintf(
      "The element `%s` in scDblFinder parameters has an invalid value.",
      broken_elem
    ))
  }
  boolean_rules <- c(
    "log_transform",
    "mean_center",
    "normalise_variance",
    "random_svd"
  )
  res <- purrr::map_lgl(boolean_rules, \(name) {
    checkmate::qtest(x[[name]], "B1")
  })
  if (!isTRUE(all(res))) {
    broken_elem <- boolean_rules[which(!res)][1]
    return(sprintf(
      "The element `%s` in scDblFinder parameters must be TRUE or FALSE.",
      broken_elem
    ))
  }
  optional_rules <- list(
    "expected_doublet_rate" = c("0", "N1(0,1]"),
    "manual_threshold" = c("0", "N1[0,)")
  )
  res <- purrr::imap_lgl(x, \(x, name) {
    if (name %in% names(optional_rules)) {
      checkmate::qtest(x, optional_rules[[name]])
    } else {
      TRUE
    }
  })
  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(sprintf(
      "The element `%s` in scDblFinder parameters is incorrect.",
      broken_elem
    ))
  }
  return(TRUE)
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
  # Checkmate extension
  res <- checkmate::checkList(x)
  if (!isTRUE(res)) {
    return(res)
  }
  res <- checkmate::checkNames(
    names(x),
    must.include = c(
      "method",
      "loess_span",
      "num_bin",
      "bin_method"
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }
  rules <- list(
    "loess_span" = "N1[0.1, 1]",
    "num_bin" = "I1"
  )
  res <- purrr::imap_lgl(x, \(x, name) {
    if (name %in% names(rules)) {
      checkmate::qtest(x, rules[[name]])
    } else {
      TRUE
    }
  })
  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(
      sprintf(
        paste(
          "The following element `%s` in single cell HVG selection is",
          "incorrect: loess_span needs to be between 0.1 and 1 and num_bin",
          "an integer."
        ),
        broken_elem
      )
    )
  }
  # test
  test_choice_rules <- list(
    method = c("vst", "meanvarbin", "dispersion"),
    bin_method = c("equal_width", "equal_freq")
  )
  test_choice_res <- purrr::imap_lgl(x, \(x, name) {
    if (name %in% names(test_choice_rules)) {
      checkmate::testChoice(x, test_choice_rules[[name]])
    } else {
      TRUE
    }
  })
  if (!isTRUE(all(test_choice_res))) {
    broken_elem <- names(test_choice_res)[which(!test_choice_res)][1]
    return(
      sprintf(
        paste0(
          "The following element `%s` in HVG params is not one of the",
          "expected choices. Please double check the documentation."
        ),
        broken_elem
      )
    )
  }

  return(TRUE)
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

#### knn -----------------------------------------------------------------------

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
  res <- checkmate::checkList(x)
  if (!isTRUE(res)) {
    return(res)
  }
  res <- checkmate::checkNames(
    names(x),
    must.include = c(
      "full_snn",
      "pruning",
      "snn_similarity",
      "ann_dist"
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }
  # KNN params
  knn_params <- x[names(x) %in% KNN_PARAM_NAMES]
  res <- checkKnnParams(knn_params)
  if (!isTRUE(res)) {
    return(res)
  }

  # Check non-kNN parameters
  rules <- list(
    "full_snn" = "B1",
    "pruning" = "N1[0, 1]"
  )

  res <- purrr::imap_lgl(x, \(val, name) {
    if (name %in% names(rules)) {
      checkmate::qtest(val, rules[[name]])
    } else {
      TRUE
    }
  })

  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(
      sprintf(
        paste(
          "The following parameter `%s` is incorrect:",
          "full_snn must be boolean, pruning must be in [0,1]."
        ),
        broken_elem
      )
    )
  }

  res <- checkmate::testChoice(x[["snn_similarity"]], c("rank", "jaccard"))
  if (!isTRUE(res)) {
    return("snn_similarity must be either 'rank' or 'jaccard'.")
  }

  return(TRUE)
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
#' @param x The `SingleCells` object to check/assert.
#' @param cell_names String. The provided cell names.
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkCellsExist <- function(x, cell_names) {
  res <- checkmate::checkClass(x, "bixverse::SingleCells")
  if (!isTRUE(res)) {
    return(res)
  }
  res <- checkmate::qtest(cell_names, "S+")
  if (!isTRUE(res)) {
    return("The cell names need be a string vector.")
  }
  all_cell_names <- get_cell_names(x)
  res <- all(cell_names %in% all_cell_names)
  if (!isTRUE(res)) {
    return(
      paste(
        "Some of the provided cell names do not exist in the object.",
        "Please check."
      )
    )
  }
  return(TRUE)
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

#' Check metacell generation parameters
#'
#' @description Checkmate extension for checking the metacell generation
#' parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkScMetacells <- function(x) {
  res <- checkmate::checkList(x)
  if (!isTRUE(res)) {
    return(res)
  }

  res <- checkmate::checkNames(
    names(x),
    must.include = c(
      "max_shared",
      "target_no_metacells",
      "max_iter"
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

  # integers
  integer_rules <- list(
    "max_shared" = "I1[1,)",
    "target_no_metacells" = "I1[1,)",
    "max_iter" = "I1[1,)"
  )

  res <- purrr::imap_lgl(x, \(val, name) {
    if (name %in% names(integer_rules)) {
      checkmate::qtest(val, integer_rules[[name]])
    } else {
      TRUE
    }
  })

  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(
      sprintf(
        paste(
          "The following element `%s` in metacell generation is incorrect:",
          "max_shared, target_no_metacells and max_iter need to be integers >= 1."
        ),
        broken_elem
      )
    )
  }

  return(TRUE)
}

#' Assert metacell generation parameters
#'
#' @description Checkmate extension for assert the metacell generation
#' parameters.
#'
#' @inheritParams checkScMetacells
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertScMetacells <- checkmate::makeAssertionFunction(checkScMetacells)

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
  res <- checkmate::checkList(x)
  if (!isTRUE(res)) {
    return(res)
  }

  res <- checkmate::checkNames(
    names(x),
    must.include = c(
      "n_sea_cells",
      "max_fw_iters",
      "convergence_epsilon",
      "max_iter",
      "min_iter",
      "greedy_threshold",
      "graph_building",
      "pruning",
      "pruning_threshold"
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

  # Check non-kNN integer parameters
  integer_rules <- list(
    "n_sea_cells" = "I1[1,)",
    "max_fw_iters" = "I1[1,)",
    "max_iter" = "I1[1,)",
    "min_iter" = "I1[1,)",
    "greedy_threshold" = "I1[1,)"
  )

  res <- purrr::imap_lgl(x, \(val, name) {
    if (name %in% names(integer_rules)) {
      checkmate::qtest(val, integer_rules[[name]])
    } else {
      TRUE
    }
  })

  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(
      sprintf(
        paste(
          "The following element `%s` in SEACells parameters is incorrect:",
          "n_sea_cells, max_fw_iters, max_iter, min_iter and greedy_threshold",
          "need to be integers >= 1."
        ),
        broken_elem
      )
    )
  }

  # Check numeric parameters
  numeric_rules <- list(
    "convergence_epsilon" = "N1",
    "pruning_threshold" = "N1"
  )

  res <- purrr::imap_lgl(x, \(val, name) {
    if (name %in% names(numeric_rules)) {
      checkmate::qtest(val, numeric_rules[[name]])
    } else {
      TRUE
    }
  })

  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(
      sprintf(
        paste(
          "The following element `%s` in SEACells parameters is incorrect:",
          "convergence_epsilon and pruning_threshold need to be numeric."
        ),
        broken_elem
      )
    )
  }

  # Check boolean parameters
  res <- checkmate::qtest(x[["pruning"]], "B1")
  if (!isTRUE(res)) {
    return("pruning needs to be a boolean.")
  }

  # Check string parameters
  res <- checkmate::qtest(x[["graph_building"]], "S1")
  if (!isTRUE(res)) {
    return("graph_building needs to be a string.")
  }

  TRUE
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
  res <- checkmate::checkList(x)
  if (!isTRUE(res)) {
    return(res)
  }
  res <- checkmate::checkNames(
    names(x),
    must.include = c(
      "walk_length",
      "graining_factor",
      "linkage_dist"
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

  # Check non-kNN integer parameters
  integer_rules <- list(
    "walk_length" = "I1[1,)"
  )

  res <- purrr::imap_lgl(x, \(val, name) {
    if (name %in% names(integer_rules)) {
      checkmate::qtest(val, integer_rules[[name]])
    } else {
      TRUE
    }
  })

  if (!isTRUE(all(res))) {
    return("walk_length needs to be an integer >= 1.")
  }

  # Check numeric parameters
  res <- checkmate::qtest(x[["graining_factor"]], "N1")
  if (!isTRUE(res)) {
    return("graining_factor needs to be numeric.")
  }

  # Check choice parameters
  res <- checkmate::testChoice(x[["linkage_dist"]], c("complete", "average"))
  if (!isTRUE(res)) {
    return("linkage_dist must be either 'complete' or 'average'.")
  }

  return(TRUE)
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
  res <- checkmate::checkList(x)
  if (!isTRUE(res)) {
    return(res)
  }

  res <- checkmate::checkNames(
    names(x),
    must.include = c(
      "neighbours_within_batch",
      "set_op_mix_ratio",
      "local_connectivity",
      "trim"
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
  res <- checkmate::qtest(x[["neighbours_within_batch"]], "I1[1,)")
  if (!isTRUE(res)) {
    return("neighbours_within_batch needs to be an integer >= 1.")
  }

  res <- checkmate::qtest(x[["trim"]], c("0", "I1[1,)"))
  if (!isTRUE(res)) {
    return("trim needs to be NULL or an integer >= 1.")
  }

  # Check numeric parameters
  res <- checkmate::qtest(x[["set_op_mix_ratio"]], "N1[0,1]")
  if (!isTRUE(res)) {
    return("set_op_mix_ratio needs to be numeric between 0 and 1.")
  }

  res <- checkmate::qtest(x[["local_connectivity"]], "N1")
  if (!isTRUE(res)) {
    return("local_connectivity needs to be numeric.")
  }

  return(TRUE)
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
  res <- checkmate::checkList(x)
  if (!isTRUE(res)) {
    return(res)
  }
  res <- checkmate::checkNames(
    names(x),
    must.include = c(
      "ndist",
      "cos_norm",
      "no_pcs",
      "random_svd"
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
  res <- checkmate::qtest(x[["no_pcs"]], "I1[1,)")
  if (!isTRUE(res)) {
    return("no_pcs needs to be an integer >= 1.")
  }
  # Check numeric parameters
  res <- checkmate::qtest(x[["ndist"]], "N1(0,)")
  if (!isTRUE(res)) {
    return("ndist needs to be a positive numeric.")
  }
  # Check logical parameters
  logical_params <- c("cos_norm", "random_svd")
  res <- purrr::map_lgl(logical_params, \(param) {
    checkmate::qtest(x[[param]], "B1")
  })
  if (!isTRUE(all(res))) {
    broken_param <- logical_params[which(!res)][1]
    return(sprintf(
      "%s needs to be logical.",
      broken_param
    ))
  }
  return(TRUE)
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

### VISION ---------------------------------------------------------------------

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
  res <- checkmate::checkList(x)
  if (!isTRUE(res)) {
    return(res)
  }
  res <- checkmate::checkNames(
    names(x),
    must.include = c(
      "n_perm",
      "n_cluster"
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

  # Check non-kNN integer parameters
  integer_rules <- list(
    "n_perm" = "I1[1,)",
    "n_cluster" = "I1[1,)"
  )

  res <- purrr::imap_lgl(x, \(val, name) {
    if (name %in% names(integer_rules)) {
      checkmate::qtest(val, integer_rules[[name]])
    } else {
      TRUE
    }
  })

  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(
      sprintf(
        paste(
          "The following element `%s` in VISION parameters is incorrect:",
          "n_perm and n_cluster need to be integers >= 1."
        ),
        broken_elem
      )
    )
  }

  return(TRUE)
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
  res <- checkmate::checkList(x)
  if (!isTRUE(res)) {
    return(res)
  }

  res <- checkmate::checkNames(
    names(x),
    must.include = c(
      "model",
      "normalise"
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

  # Check choice parameters
  res <- checkmate::testChoice(x[["model"]], c("danb", "bernoulli", "normal"))
  if (!isTRUE(res)) {
    return("model must be one of: danb, bernoulli, normal.")
  }

  # Check boolean parameters
  res <- checkmate::qtest(x[["normalise"]], "B1")
  if (!isTRUE(res)) {
    return("normalise needs to be a boolean.")
  }

  return(TRUE)
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
  res <- checkmate::checkList(x)
  if (!isTRUE(res)) {
    return(res)
  }

  res <- checkmate::checkNames(
    names(x),
    must.include = c(
      "prop",
      "k_refine",
      "refinement_strategy",
      "index_type"
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

  # Check non-kNN integer parameters
  res <- checkmate::qtest(x[["k_refine"]], "I1[1,)")
  if (!isTRUE(res)) {
    return("k_refine must be an integer >= 1.")
  }

  # Check numeric parameters
  res <- checkmate::qtest(x[["prop"]], "N1(0,1)")
  if (!isTRUE(res)) {
    return("prop must be in (0,1).")
  }

  # Check choice parameters
  test_choice_rules <- list(
    refinement_strategy = c("approximate", "bruteforce", "index"),
    index_type = c("annoy", "hnsw")
  )

  test_choice_res <- purrr::imap_lgl(x, \(val, name) {
    if (name %in% names(test_choice_rules)) {
      checkmate::testChoice(val, test_choice_rules[[name]])
    } else {
      TRUE
    }
  })

  if (!isTRUE(all(test_choice_res))) {
    broken_elem <- names(test_choice_res)[which(!test_choice_res)][1]
    return(
      sprintf(
        paste(
          "The following element `%s` in MiloR parameters is not one of the",
          "expected choices. Please check the documentation."
        ),
        broken_elem
      )
    )
  }

  return(TRUE)
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

  # Integer rules
  integer_rules <- list(
    "k" = c("I1[1,)", "0"),
    "max_iter_kmeans" = "I1[1,)",
    "max_iter_harmony" = "I1[1,)",
    "window_size" = "I1[1,)"
  )

  res <- purrr::imap_lgl(x, \(x, name) {
    if (name %in% names(integer_rules)) {
      checkmate::qtest(x, integer_rules[[name]])
    } else {
      TRUE
    }
  })

  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(
      sprintf(
        paste(
          "The following element `%s` in Harmony parameters is incorrect:",
          "max_iter_kmeans, max_iter_harmony,",
          "and window_size must be integers >= 1. k must be NULL or an integer."
        ),
        broken_elem
      )
    )
  }

  # Numeric vector rules (can be length 1 or longer)
  vector_rules <- list(
    "sigma" = "N+[0,)",
    "theta" = "N+[0,)",
    "lambda" = "N+[0,)"
  )

  res <- purrr::imap_lgl(x, \(x, name) {
    if (name %in% names(vector_rules)) {
      checkmate::qtest(x, vector_rules[[name]])
    } else {
      TRUE
    }
  })

  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(
      sprintf(
        paste(
          "The following element `%s` in Harmony parameters is incorrect:",
          "sigma, theta, and lambda must be numeric vectors",
          "with non-negative values."
        ),
        broken_elem
      )
    )
  }

  # Scalar numeric rules
  scalar_rules <- list(
    "block_size" = "N1(0,1]",
    "epsilon_kmeans" = "N1(0,)",
    "epsilon_harmony" = "N1(0,)"
  )

  res <- purrr::imap_lgl(x, \(x, name) {
    if (name %in% names(scalar_rules)) {
      checkmate::qtest(x, scalar_rules[[name]])
    } else {
      TRUE
    }
  })

  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(
      sprintf(
        paste(
          "The following element `%s` in Harmony parameters is incorrect:",
          "block_size must be in (0,1]; epsilon_kmeans",
          "and epsilon_harmony must be > 0."
        ),
        broken_elem
      )
    )
  }

  return(TRUE)
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
  res <- checkmate::checkList(x)
  if (!isTRUE(res)) {
    return(res)
  }

  # Fields shared across all learner types
  required_names <- c(
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
  res <- checkmate::checkNames(names(x), must.include = required_names)
  if (!isTRUE(res)) {
    return(res)
  }

  res <- checkmate::checkChoice(
    x$learner_type,
    c("randomforest", "extratrees", "grnboost2")
  )
  if (!isTRUE(res)) {
    return(paste("learner_type:", res))
  }

  res <- checkmate::checkChoice(
    x$gene_batch_strategy,
    c("random", "correlated")
  )
  if (!isTRUE(res)) {
    return(paste("gene_batch_strategy:", res))
  }

  # Integer validation shared across all types
  integer_rules <- list(
    min_counts = "I1[1,)",
    n_pcs = "I1[1,)",
    n_subsample = "I1[1,)",
    min_samples_leaf = "I1[1,)",
    n_features_split = "I1[0,)",
    max_depth = "I1[1,)"
  )
  res <- purrr::imap_lgl(x, \(val, name) {
    if (name %in% names(integer_rules)) {
      checkmate::qtest(val, integer_rules[[name]])
    } else {
      TRUE
    }
  })
  if (!all(res)) {
    return(sprintf(
      "Element `%s` in SCENIC parameters failed integer validation.",
      names(res)[!res][1]
    ))
  }

  # Numeric validation shared across all types
  scalar_numeric_rules <- list(
    min_cells = "N1(0,1]"
  )
  res <- purrr::imap_lgl(x, \(val, name) {
    if (name %in% names(scalar_numeric_rules)) {
      checkmate::qtest(val, scalar_numeric_rules[[name]])
    } else {
      TRUE
    }
  })
  if (!all(res)) {
    return(sprintf(
      "Element `%s` in SCENIC parameters failed numeric validation.",
      names(res)[!res][1]
    ))
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
    if (!is.null(x$subsample_frac)) {
      if (!checkmate::qtest(x$subsample_frac, "N1(0,1]")) {
        return("subsample_frac must be a numeric in (0, 1] or NULL.")
      }
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
    if (!is.null(x$subsample_frac)) {
      if (!checkmate::qtest(x$subsample_frac, "N1(0,1]")) {
        return("subsample_frac must be a numeric in (0, 1] or NULL.")
      }
    }
  }

  if (x$learner_type == "grnboost2") {
    if (
      is.null(x$n_trees_max) ||
        !checkmate::qtest(x$n_trees_max, "I1[1,)")
    ) {
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
