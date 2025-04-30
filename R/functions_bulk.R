# helpers ----------------------------------------------------------------------

## dge helpers -----------------------------------------------------------------

#' Fixes contrast names for DGEs
#'
#' @param x Vector of strings or factors.
#'
#' @returns Vector with fixed naming based on R conventions.
#'
#' @export
fix_contrast_names <- function(x) {
  checkmate::qassert(x, c("S+", "F+", "N+"))
  if (checkmate::qtest(x, c("S+", "F+"))) {
    res <- as.factor(gsub("\\.", "_", make.names(gsub("[[:punct:]]", "", x))))
  } else {
    res <- x
  }
  return(res)
}


#' Create all limma contrasts
#'
#' @param limma_fit The fitted limma model, i.e., output of [limma::lmFit()].
#'
#' @returns The Limma contrasts for further usage.
#'
#' @export
all_limma_contrasts <- function(limma_fit) {
  checkmate::assertClass(limma_fit, "MArrayLM")
  coefs_fit <- colnames(coef(limma_fit))
  cb <- combn(x = coefs_fit, m = 2, FUN = function(x) {
    paste0(x[[1]], "-", x[[2]])
  })
  contrasts <- limma::makeContrasts(contrasts = cb, levels = coefs_fit)
  return(contrasts)
}

# dge functions ----------------------------------------------------------------

## traditional dge -------------------------------------------------------------

#' Wrapper for a Limma Voom analysis
#'
#' @description
#' Wrapper function to run Limma Voom workflows.
#'
#' @param meta_data data.table. The meta information about the experiment in
#' which the contrast info (and potential co-variates) can be found.
#' @param main_contrast String. Which column contains the main groups you want
#' to test differential gene expression with the Limma-Voom workflow for.
#' @param normalized_counts matrix, normalized count matrix.
#' @param co_variates String or NULL. Optional co-variates you wish to consider
#' during model fitting.
#' @param ... Additional parameters to forward to [edgeR::filterByExpr()].
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @returns A data.table with all the DGE results from [limma::topTable()] for
#' the identified contrast pairs.
#'
#' @export
#'
#' @import data.table
#' @importFrom magrittr `%>%`
run_limma_voom <- function(
  meta_data,
  main_contrast,
  normalized_counts,
  co_variates = NULL,
  ...,
  .verbose = TRUE
) {
  variables <- c(main_contrast, co_variates)
  # Checks
  checkmate::assertDataFrame(meta_data)
  checkmate::assertNames(
    names(S7::prop(object, "outputs")),
    must.include = "normalized_counts"
  )
  checkmate::qassert(co_variates, c("S+", "0"))
  checkmate::assertNames(
    names(meta_data),
    must.include = variables
  )
  checkmate::qassert(.verbose, "B1")

  # Fix any names
  meta_data[,
    (variables) := lapply(.SD, fix_contrast_names),
    .SDcols = variables
  ]
  if (.verbose)
    message(paste(
      "Fixing any naming issues for the selected main contrast",
      "and any co-variates."
    ))

  model_formula <- sprintf(
    "~ 0 + %s",
    paste(variables, collapse = " + ")
  )

  model_matrix <- model.matrix(as.formula(model_formula), data = meta_data)
  colnames(model_matrix) <- gsub(main_contrast, "", colnames(model_matrix))

  limma_fit <- limma::lmFit(
    normalized_counts[, meta_data$sample_id],
    model_matrix
  )

  contrasts <- all_limma_contrasts(limma_fit)

  final_fit <- limma::contrasts.fit(limma_fit, contrasts)
  final_fit <- limma::eBayes(final_fit)

  tested_contrasts <- attributes(contrasts)$dimnames$Contrasts
  if (.verbose)
    message(paste(
      "Tested these contrasts:",
      paste(tested_contrasts, collapse = ", ")
    ))

  all_dge_res <- purrr::map(tested_contrasts, \(coef) {
    top.table <- as.data.table(
      limma::topTable(
        fit = final_fit,
        coef = coef,
        sort.by = "P",
        n = Inf,
        confint = TRUE
      ),
      keep.rownames = "gene_id"
    ) %>%
      .[, contrast := gsub("-", "_vs_", coef)]
  }) %>%
    rbindlist()

  return(all_dge_res)
}

## effect size calculations ----------------------------------------------------

#' Calculate the effect
#'
#' @param meta_data data.table. The meta information about the experiment in
#' which the contrast info can be found.
#' @param main_contrast String. Which column contains the main groups you want
#' to calculate the Hedge's G effect for. Every permutation of the groups
#' will be tested.
#' @param normalized_counts normalized count matrix.
#' @param ... Additional parameters to forward to [edgeR::filterByExpr()].
#' @param .verbose Controls verbosity of the function.
#'
#' @returns A data.table with the effect sizes and standard errors based on the
#' Hedge's G effect size for all found permutations of the groups.
#'
#' @export
#'
#' @import data.table
#' @importFrom magrittr `%>%`
hedges_g_dge_list <- function(
  meta_data,
  main_contrast,
  normalized_counts,
  ...,
  .verbose = TRUE
) {
  # Checks
  checkmate::assertDataFrame(meta_data)
  checkmate::qassert(main_contrast, "S1")
  checkmate::assertClass(normalized_counts, "matrix")
  checkmate::assertNames(
    names(meta_data),
    must.include = main_contrast
  )
  # Function
  groups <- as.character(unique(meta_data[[main_contrast]]))
  combinations_to_test <- combn(
    x = groups,
    m = 2,
    FUN = function(x) {
      c(x[[1]], x[[2]])
    },
    simplify = FALSE
  )

  res <- purrr::map(combinations_to_test, \(combination) {
    grpA <- meta_data[
      eval(parse(text = paste0(main_contrast, " == '", combination[[1]], "'"))),
      sample_id
    ]
    grpB <- meta_data[
      eval(parse(text = paste0(main_contrast, " == '", combination[[2]], "'"))),
      sample_id
    ]

    mat_a <- t(normalized_counts[, grpA])
    mat_b <- t(normalized_counts[, grpB])

    hedges_g_effect <- calculate_effect_size(
      mat_a = mat_a,
      mat_b = mat_b,
      .verbose = .verbose
    ) %>%
      data.table::setDT() %>%
      .[, `:=`(
        gene_id = colnames(mat_a),
        combination = paste(combination[[1]], combination[[2]], sep = "_vs_")
      )]

    hedges_g_effect
  }) %>%
    data.table::rbindlist()

  return(res)
}
