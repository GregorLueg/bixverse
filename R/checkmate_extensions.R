# checks -----------------------------------------------------------------------

## others ----------------------------------------------------------------------

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

#' Check ssGSEA parameters
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

## dgrdl -----------------------------------------------------------------------

#' Check DGRDL parameters
#'
#' @description Checkmate extension for checking dual graph regularised
#' dictionary learning parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
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

## snf -------------------------------------------------------------------------

#' Check SNF parameters
#'
#' @description Checkmate extension for checking the SNF parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
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

## cistarget -------------------------------------------------------------------

#' Check CisTarget parameters
#'
#' @description Checkmate extension for checking CisTarget parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
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

## single cell -----------------------------------------------------------------

### synthetic data -------------------------------------------------------------
#' Check synthetic data parameters
#'
#' @description Checkmate extension for checking the synthetic data
#' parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
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
      "batch_effect_strength"
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  rules <- list(
    "n_cells" = "I1",
    "n_genes" = "I1",
    "n_batches" = "I1"
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
          "n_cells, n_genes and n_batches need to be integers."
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

  return(TRUE)
}

### io -------------------------------------------------------------------------

#' Check SC MTX load parameters
#'
#' @description Checkmate extension for checking MTX loading parameters
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
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

### qc -------------------------------------------------------------------------

#' Check SC minimum QC parameters
#'
#' @description Checkmate extension for checking the minimum QC parameters
#' in single cell.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
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

### scrublet -------------------------------------------------------------------

#' Check Scrublet parameters
#'
#' @description Checkmate extension for checking the Scrublet parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
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
      "random_svd",
      "k",
      "knn_method",
      "dist_metric",
      "search_budget",
      "n_trees",
      "nn_max_iter",
      "rho",
      "delta"
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  # Integer rules
  integer_rules <- list(
    "no_pcs" = "I1[1,)",
    "k" = "I1[0,)",
    "search_budget" = "I1[1,)",
    "n_trees" = "I1[1,)",
    "n_bins" = "I1[10,)",
    "nn_max_iter" = "I1[1,)"
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
          "no_pcs must be >= 1; k must be >= 0;",
          "search_budget, n_trees and nn_max_iter must be >= 1;",
          "n_bins must be >= 10."
        ),
        broken_elem
      )
    )
  }

  # Numeric rules
  numeric_rules <- list(
    "min_gene_var_pctl" = "N1[0,1]",
    "loess_span" = "N1(0,)",
    "sim_doublet_ratio" = "N1(0,)",
    "expected_doublet_rate" = "N1[0,1]",
    "stdev_doublet_rate" = "N1[0,1]",
    "rho" = "N1(0,)",
    "delta" = "N1[0,1]",
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
          "min_gene_var_pctl, expected_doublet_rate, stdev_doublet_rate",
          "and delta must be in [0,1];",
          "loess_span, rho sim_doublet_ratio must be > 0;",
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
  test_choice_rules <- list(
    hvg_method = c("vst", "mvb", "dispersion"),
    knn_method = c("annoy", "hnsw", "nndescent"),
    dist_metric = c("euclidean", "cosine")
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
          "The following element `%s` in the Scrublet parameters is not one of",
          " the expected choices. Please double check the documentation."
        ),
        broken_elem
      )
    )
  }

  return(TRUE)
}

### boost ----------------------------------------------------------------------

#' Check Boost parameters
#'
#' @description Checkmate extension for checking the Boost parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
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
      "voter_thresh",
      "k",
      "knn_method",
      "dist_metric",
      "search_budget",
      "n_trees",
      "nn_max_iter",
      "rho",
      "delta"
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  # Integer rules
  integer_rules <- list(
    "no_pcs" = "I1[1,)",
    "n_iters" = "I1[1,)",
    "k" = "I1[0,)",
    "search_budget" = "I1[1,)",
    "n_trees" = "I1[1,)",
    "nn_max_iter" = "I1[1,)"
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
          "no_pcs, n_iters, search_budget, nn_max_iter & n_trees must be >= 1;",
          "k must be >= 0."
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
    "rho" = "N1(0,)",
    "delta" = "N1[0,1]",
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
          "min_gene_var_pctl, boost_rate, delta and voter_thresh must be in",
          "[0,1]; loess_span, resolution, rho and p_thresh must be > 0;",
          "target_size must be > 1"
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
          "clip_max and target_size must be NULL or positive numeric values."
        ),
        broken_elem
      )
    )
  }

  # Choice rules
  test_choice_rules <- list(
    hvg_method = c("vst", "mvb", "dispersion"),
    knn_method = c("annoy", "hnsw", "nndescent"),
    dist_metric = c("euclidean", "cosine")
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
          "The following element `%s` in the Boost parameters is not one of",
          " the expected choices. Please double check the documentation."
        ),
        broken_elem
      )
    )
  }

  return(TRUE)
}

### hvg ------------------------------------------------------------------------

#' Check HVG selection parameters
#'
#' @description Checkmate extension for checking the HVG parameters for single
#' cell.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
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

### knn ------------------------------------------------------------------------

#' Check neighbour generation parameters
#'
#' @description Checkmate extension for checking the neighbour generation
#' parameters for single cell.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
checkScNeighbours <- function(x) {
  res <- checkmate::checkList(x)
  if (!isTRUE(res)) {
    return(res)
  }
  res <- checkmate::checkNames(
    names(x),
    must.include = c(
      "k",
      "n_trees",
      "search_budget",
      "max_iter",
      "rho",
      "delta",
      "knn_algorithm",
      "full_snn",
      "pruning",
      "snn_similarity",
      "ann_dist"
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }
  rules <- list(
    "k" = "I1",
    "n_trees" = "I1",
    "search_budget" = "I1",
    "max_iter" = "I1",
    "rho" = "N1",
    "delta" = "N1",
    "full_snn" = "B1",
    "pruning" = "N1[0, 1]"
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
          "The following element `%s` in single cell KNN generation is",
          "incorrect: k, n_trees, search_budget, and max_iter need to be integers.",
          "rho and delta need to be numeric. full_snn needs to be boolean and",
          "pruning a number between 0 and 1."
        ),
        broken_elem
      )
    )
  }
  test_choice_rules <- list(
    knn_algorithm = c("annoy", "hnsw", "nndescent"),
    snn_similarity = c("rank", "jaccard"),
    ann_dist = c("cosine", "euclidean")
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
          "The following element `%s` in the KNN generation is not one of ",
          "the expected choices. Please double check the documentation."
        ),
        broken_elem
      )
    )
  }
  return(TRUE)
}

### cells in object ------------------------------------------------------------

#' Check that the cell name exists in the object
#'
#' @description Checkmate extension for checking if the prodivided cell names
#' exist in the object.
#'
#' @param x The `single_cell_exp` object to check/assert.
#' @param cell_names String. The provided cell names.
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
checkCellsExist <- function(x, cell_names) {
  res <- checkmate::checkClass(x, "bixverse::single_cell_exp")
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

### meta cells -----------------------------------------------------------------

#' Check metacell generation parameters
#'
#' @description Checkmate extension for checking the metacell generation
#' parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
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
      "max_iter",
      "k",
      "knn_method",
      "n_trees",
      "search_budget",
      "ann_dist"
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  rules <- list(
    "max_shared" = "I1",
    "target_no_metacells" = "I1",
    "max_iter" = "I1",
    "k" = "I1",
    "n_trees" = "I1",
    "search_budget" = "I1"
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
          "The following element `%s` in metacell generation is incorrect:",
          "max_shared, target_no_metacells, max_iter, k, n_trees and",
          "search_budget need to be integers."
        ),
        broken_elem
      )
    )
  }

  test_choice_rules <- list(
    knn_method = c("annoy", "hnsw"),
    ann_dist = c("cosine", "euclidean")
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
          "The following element `%s` in the metacell generation is not one of",
          " the expected choices. Please double check the documentation."
        ),
        broken_elem
      )
    )
  }

  return(TRUE)
}

### bbknn ----------------------------------------------------------------------

#' Check BBKNN parameters
#'
#' @description Checkmate extension for checking the BBKNN parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
checkScBbknn <- function(x) {
  res <- checkmate::checkList(x)
  if (!isTRUE(res)) {
    return(res)
  }

  res <- checkmate::checkNames(
    names(x),
    must.include = c(
      "neighbours_within_batch",
      "knn_method",
      "ann_dist",
      "set_op_mix_ratio",
      "local_connectivity",
      "annoy_n_trees",
      "search_budget",
      "trim"
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  integer_rules <- list(
    "neighbours_within_batch" = "I1",
    "annoy_n_trees" = "I1",
    "search_budget" = "I1",
    "trim" = c("0", "I1")
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
          "The following element `%s` in BBKNN parameters is incorrect:",
          "neighbours_within_batch, annoy_n_trees, and search_budget need to be integers;",
          "trim needs to be NULL or an integer."
        ),
        broken_elem
      )
    )
  }

  numeric_rules <- list(
    "set_op_mix_ratio" = "N1[0, 1]",
    "local_connectivity" = "N1"
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
          "The following element `%s` in BBKNN parameters is incorrect:",
          "set_op_mix_ratio and local_connectivity need to be numeric."
        ),
        broken_elem
      )
    )
  }

  test_choice_rules <- list(
    knn_method = c("annoy", "hnsw"),
    ann_dist = c("cosine", "euclidean")
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
          "The following element `%s` in the BBKNN parameters is not one of",
          " the expected choices. Please double check the documentation."
        ),
        broken_elem
      )
    )
  }

  return(TRUE)
}

### fastMNN --------------------------------------------------------------------

#' Check fastMNN parameters
#'
#' @description Checkmate extension for checking the fastMNN parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
checkScFastmnn <- function(x) {
  res <- checkmate::checkList(x)
  if (!isTRUE(res)) {
    return(res)
  }

  res <- checkmate::checkNames(
    names(x),
    must.include = c(
      "k",
      "sigma",
      "knn_method",
      "dist_metric",
      "annoy_n_trees",
      "annoy_search_budget",
      "cos_norm",
      "var_adj",
      "no_pcs",
      "random_svd"
    )
  )
  if (!isTRUE(res)) {
    return(res)
  }

  integer_rules <- list(
    "k" = "I1",
    "annoy_n_trees" = "I1",
    "annoy_search_budget" = "I1",
    "no_pcs" = "I1"
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
      paste(
        "The following element `%s` in fastMNN parameters is incorrect: k,",
        "annoy_n_trees, annoy_search_budget, and no_pcs need to be integers."
      ),
      broken_elem
    ))
  }

  numeric_rules <- list("sigma" = "N1")
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
      paste(
        "The following element `%s` in fastMNN parameters is incorrect:",
        "sigma needs to be numeric."
      ),
      broken_elem
    ))
  }

  logical_rules <- list(
    "cos_norm" = "B1",
    "var_adj" = "B1",
    "random_svd" = "B1"
  )
  res <- purrr::imap_lgl(x, \(x, name) {
    if (name %in% names(logical_rules)) {
      checkmate::qtest(x, logical_rules[[name]])
    } else {
      TRUE
    }
  })
  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(sprintf(
      paste(
        "The following element `%s` in fastMNN parameters is incorrect:",
        "cos_norm, var_adj, and random_svd need to be logical."
      ),
      broken_elem
    ))
  }

  test_choice_rules <- list(
    knn_method = c("annoy", "hnsw"),
    dist_metric = c("cosine", "euclidean")
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
    return(sprintf(
      paste(
        "The following element `%s` in fastMNN parameters is not one of the",
        "expected choices. Please check the documentation."
      ),
      broken_elem
    ))
  }

  TRUE
}

# asserts ----------------------------------------------------------------------

## other -----------------------------------------------------------------------

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
assertFileExists <- checkmate::makeAssertionFunction(checkFilesExist)

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

## ssgsea ----------------------------------------------------------------------

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
assertSingleSampleGSEAparams <- checkmate::makeAssertionFunction(
  checkSingleSampleGSEAparams
)

## gsva ------------------------------------------------------------------------

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

## DGRDL ------------------------------------------------------------------------

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
assertDGRDLparams <- checkmate::makeAssertionFunction(checkDGRDLparams)

## SNF -------------------------------------------------------------------------

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
assertCistargetParams <- checkmate::makeAssertionFunction(checkCistargetParams)

## CisTarget -------------------------------------------------------------------

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
assertSNFParams <- checkmate::makeAssertionFunction(checkSNFParams)

## single cell -----------------------------------------------------------------

### synthetic data -------------------------------------------------------------

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
assertScSyntheticData <- checkmate::makeAssertionFunction(checkScSyntheticData)

### io -------------------------------------------------------------------------

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
assertScMtxIO <- checkmate::makeAssertionFunction(checkScMtxIO)

### qc -------------------------------------------------------------------------

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
assertScMinQC <- checkmate::makeAssertionFunction(checkScMinQC)

### scrublet -------------------------------------------------------------------

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
assertScScrublet <- checkmate::makeAssertionFunction(checkScScrublet)

### boost ----------------------------------------------------------------------

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
assertScBoost <- checkmate::makeAssertionFunction(checkScBoost)

### hvg ------------------------------------------------------------------------

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
assertScHvg <- checkmate::makeAssertionFunction(checkScHvg)

### knn ------------------------------------------------------------------------

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
assertScNeighbours <- checkmate::makeAssertionFunction(checkScNeighbours)

### cell exists ----------------------------------------------------------------

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
assertCellsExist <- checkmate::makeAssertionFunction(checkCellsExist)

### meta cells -----------------------------------------------------------------

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
assertScMetacells <- checkmate::makeAssertionFunction(checkScMetacells)

### bbknn ----------------------------------------------------------------------

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
assertScBbknn <- checkmate::makeAssertionFunction(checkScBbknn)

### fastMNN --------------------------------------------------------------------

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
assertScFastmnn <- checkmate::makeAssertionFunction(checkScFastmnn)

# tests ------------------------------------------------------------------------

## SNF -------------------------------------------------------------------------

#' @rdname checkSNFParams
#'
#' @export
testSNFParams <- checkmate::makeTestFunction(checkSNFParams)
