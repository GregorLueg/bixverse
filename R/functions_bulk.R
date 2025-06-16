# helpers ----------------------------------------------------------------------

## dge helpers -----------------------------------------------------------------

#' Fixes contrast names for DGEs
#'
#' @param x Vector of strings or factors.
#'
#' @returns Vector with fixed naming based on R conventions.
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
#' @param contrast_grps String vector. The contrast groups of interest. If NULL
#' all co-variate comparisons will be returned.
#'
#' @returns The Limma contrasts for further usage.
all_limma_contrasts <- function(limma_fit, contrast_grps) {
  # Checks
  checkmate::assertClass(limma_fit, "MArrayLM")
  checkmate::qassert(contrast_grps, c("S+", "F+", "0"))
  coefs_fit <- colnames(coef(limma_fit))
  if (!is.null(contrast_grps)) {
    coefs_fit <- coefs_fit[coefs_fit %in% contrast_grps]
  }
  contrast_combs <- c(combn(x = coefs_fit, m = 2, FUN = function(x) {
    paste0(x[[1]], "-", x[[2]])
  }))

  all_contrasts <- vector(mode = "list", length = length(contrast_combs))

  # makeContrasts is WILD...
  for (i in seq_along(contrast_combs)) {
    value <- contrast_combs[[i]]
    contrast_i <- limma::makeContrasts(
      contrasts = value,
      levels = colnames(coef(limma_fit))
    )
    all_contrasts[[i]] <- contrast_i
  }

  return(all_contrasts)
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
#' @param dge_list DGEList, see [edgeR::DGEList()].
#' @param co_variates String or NULL. Optional co-variates you wish to consider
#' during model fitting.
#' @param ... Additional parameters to forward to [limma::eBayes()] or
#' [limma::voom()].
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
  dge_list,
  co_variates = NULL,
  ...,
  .verbose = TRUE
) {
  variables <- c(main_contrast, co_variates)
  # Checks
  checkmate::assertDataFrame(meta_data)
  checkmate::qassert(main_contrast, "S1")
  checkmate::assertClass(dge_list, "DGEList")
  checkmate::qassert(co_variates, c("S+", "0"))
  checkmate::assertNames(
    names(meta_data),
    must.include = variables
  )
  checkmate::qassert(.verbose, "B1")

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

  voom_obj <- limma::voom(
    counts = dge_list,
    design = model_matrix,
    normalize.method = "quantile",
    plot = TRUE,
    ...
  )

  limma_fit <- limma::lmFit(voom_obj, model_matrix)
  contrast_grps <- unique(meta_data[[main_contrast]])

  limma_contrasts <- all_limma_contrasts(
    limma_fit = limma_fit,
    contrast_grps = contrast_grps
  )

  all_dge_res <- purrr::map(limma_contrasts, \(contrast_obj) {
    final_fit <- limma::contrasts.fit(limma_fit, contrast_obj)
    final_fit <- limma::eBayes(final_fit, ...)

    coef_name <- colnames(coef(final_fit))

    top.table <- as.data.table(
      limma::topTable(
        fit = final_fit,
        sort.by = "P",
        n = Inf,
        confint = TRUE
      ),
      keep.rownames = "gene_id"
    ) %>%
      .[, contrast := gsub("-", "_vs_", coef_name)]
  }) %>%
    data.table::rbindlist()

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
#' @param normalised_counts Numeric Matrix. The normalized count matrix.
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @returns A data.table with the effect sizes and standard errors based on the
#' Hedge's G effect size for all found permutations of the groups.
#'
#' @export
#'
#' @import data.table
#' @importFrom magrittr `%>%`
hedges_g_dge <- function(
  meta_data,
  main_contrast,
  normalised_counts,
  .verbose = TRUE
) {
  # Checks
  checkmate::assertDataFrame(meta_data)
  checkmate::qassert(main_contrast, "S1")
  checkmate::assertClass(normalised_counts, "matrix")
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

    mat_a <- t(normalised_counts[, grpA])
    mat_b <- t(normalised_counts[, grpB])

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
