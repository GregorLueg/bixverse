# checks -----------------------------------------------------------------------

## correlation params ----------------------------------------------------------

#' Check correlation graph parameters
#'
#' @description Checkmate extension for checking the graph parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
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

#' Check resolution graph parameters
#'
#' @description Checkmate extension for checking the resolution parameters for
#' community detection with Leiden.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
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

## ica params ------------------------------------------------------------------

#' Check ICA parameters
#'
#' @description Checkmate extension for checking the ICA parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
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


#' Check ICA no of component parameters
#'
#' @description Checkmate extension for checking the ICA number of component
#' parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
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

#' Check ICA randomisation parameters
#'
#' @description Checkmate extension for checking the ICA randomisation parameters
#' for a version of stabilised ICA.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
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

## community detections --------------------------------------------------------

#' Check community detection parameters
#'
#' @description Checkmate extension for checking the community detection
#' parameters for identifying genetically privileged communities.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
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

## gsea ------------------------------------------------------------------------

#' Check GSEA parameters
#'
#' @description Checkmate extension for checking the gene set enrichment
#' analysis (GSEA) parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
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

## gsva ------------------------------------------------------------------------

#' Check GSVA parameters
#'
#' @description Checkmate extension for checking the gene set variation analysis
#' (GSVA) parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
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

## ssgsea ----------------------------------------------------------------------

#' Check GSVA parameters
#'
#' @description Checkmate extension for checking single sample gene set
#' enrichment analysis parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
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
          "The following element `%s` in GSVA params does not conform to the",
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

## coremo ----------------------------------------------------------------------

#' Check CoReMo parameters
#'
#' @description Checkmate extension for checking the CoReMo parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
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

# asserts ----------------------------------------------------------------------

## correlation params ----------------------------------------------------------

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
assertCorGraphParams <- checkmate::makeAssertionFunction(checkCorGraphParams)


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
assertGraphResParams <- checkmate::makeAssertionFunction(checkGraphResParams)

## ica params ------------------------------------------------------------------

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
assertIcaParams <- checkmate::makeAssertionFunction(checkIcaParams)

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
assertIcaNcomps <- checkmate::makeAssertionFunction(checkIcaNcomps)

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
assertIcaIterParams <- checkmate::makeAssertionFunction(checkIcaIterParams)

## community detections --------------------------------------------------------

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
assertCommunityParams <- checkmate::makeAssertionFunction(checkCommunityParams)

## gsea ------------------------------------------------------------------------

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
assertGSEAParams <- checkmate::makeAssertionFunction(checkGSEAParams)

## gsva ------------------------------------------------------------------------

#' Assert GSVA parameter
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
assertSingleSampleGSEAparams <- checkmate::makeAssertionFunction(
  checkSingleSampleGSEAparams
)

## ssgsea ----------------------------------------------------------------------

#' Assert ssGSEA parameter
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
assertGSVAParams <- checkmate::makeAssertionFunction(checkGSVAParams)

## coremo ----------------------------------------------------------------------

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
assertCoReMoParams <- checkmate::makeAssertionFunction(checkCoReMoParams)
