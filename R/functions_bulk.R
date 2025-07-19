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
    res <- as.factor(gsub(
      "\\.",
      "_",
      make.names(gsub("[[:punct:]&&[^_]]", "", x))
    ))
  } else {
    res <- x
  }
  return(res)
}

#' Create all limma contrasts (based on combination of everything)
#'
#' @param limma_fit The fitted limma model, i.e., output of [limma::lmFit()].
#' @param contrast_grps String vector. The contrast groups of interest. If NULL
#' all co-variate comparisons will be returned.
#'
#' @returns The Limma contrasts for further usage.
all_limma_contrasts <- function(limma_fit, contrast_grps) {
  # Globals
  coef <- combn <- NULL

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

#' Create all limma contrasts from a provided string
#'
#' @param limma_fit The fitted limma model, i.e., output of [limma::lmFit()].
#' @param contrast_list String vector. The strings need to have the form of
#' `"contrast1-contrast2"`.
#'
#' @returns The Limma contrasts for further usage.
limma_contrasts <- function(limma_fit, contrast_list) {
  # Checks
  checkmate::assertClass(limma_fit, "MArrayLM")
  checkmate::qassert(contrast_list, c("S+", "F+", "0"))

  coefs_fit <- colnames(coef(limma_fit))
  all_contrasts <- vector(mode = "list", length = length(contrast_list))

  # makeContrasts is WILD...
  for (i in seq_along(contrast_list)) {
    value <- contrast_list[[i]]
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
#' @param contrast_list String vector or NULL. Optional string vector of
#' contrast formatted as `"contrast1-contrast2"`. Default NULL will create all
#' contrasts automatically.
#' @param dge_list DGEList, see [edgeR::DGEList()]. If `NULL` it will be
#' regenerated.
#' @param normalised_counts Matrix. Normalised count matrix.
#' @param co_variates String or NULL. Optional co-variates you wish to consider
#' during model fitting.
#' @param parallel Boolean. Shall parallelisation be used across the contrasts.
#' @param no_cores Optional integer. Number of cores to use for parallelisation.
#' If `NULL` will default to half the detected cores.
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
  contrast_list = NULL,
  dge_list = NULL,
  normalised_counts = NULL,
  co_variates = NULL,
  parallel = TRUE,
  no_cores = NULL,
  ...,
  .verbose = TRUE
) {
  variables <- c(main_contrast, co_variates)
  # Checks
  checkmate::assertDataFrame(meta_data)
  checkmate::qassert(main_contrast, "S1")
  checkmate::assert(
    checkmate::checkClass(dge_list, "DGEList"),
    checkmate::checkNull(dge_list),
    combine = "or"
  )
  checkmate::qassert(co_variates, c("S+", "0"))
  checkmate::assertNames(
    names(meta_data),
    must.include = variables
  )
  checkmate::qassert(.verbose, "B1")
  checkmate::qassert(contrast_list, c("S+", "0"))
  checkmate::assert(all(grepl("-", contrast_list)))
  checkmate::assertMatrix(normalised_counts)
  checkmate::qassert(parallel, "B1")
  checkmate::qassert(no_cores, c("I1", "0"))

  meta_data[,
    (variables) := lapply(.SD, fix_contrast_names),
    .SDcols = variables
  ]
  if (.verbose) {
    message(paste(
      "Fixing any naming issues for the selected main contrast",
      "and any co-variates."
    ))
  }

  model_formula <- sprintf(
    "~ 0 + %s",
    paste(variables, collapse = " + ")
  )

  model_matrix <- model.matrix(as.formula(model_formula), data = meta_data)
  colnames(model_matrix) <- gsub(main_contrast, "", colnames(model_matrix))

  if (is.null(dge_list)) {
    counts = normalised_counts
    limma_fit <- limma::lmFit(counts, model_matrix)
  } else {
    counts = dge_list
    voom_obj <- limma::voom(
      counts = counts,
      design = model_matrix,
      normalize.method = "quantile",
      plot = TRUE,
      ...
    )
    limma_fit <- limma::lmFit(voom_obj, model_matrix)
  }

  contrast_grps <- unique(meta_data[[main_contrast]])

  if (is.null(contrast_list)) {
    limma_contrasts <- all_limma_contrasts(
      limma_fit = limma_fit,
      contrast_grps = contrast_grps
    )
  } else {
    limma_contrasts <- limma_contrasts(
      limma_fit = limma_fit,
      contrast_list = contrast_list
    )
  }

  # parallelisation stuff
  if (!parallel) {
    mirai::daemons(1)
  } else {
    sessions <- ifelse(is.null(no_cores), get_cores(), no_cores)
    mirai::daemons(sessions)
  }

  all_dge_res <- mirai::mirai_map(
    limma_contrasts,
    \(contrast_obj, limma_fit) {
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
    },
    .args = list(limma_fit = limma_fit)
  )[] %>%
    rbindlist()

  mirai::daemons(0)

  return(all_dge_res)
}

## effect size calculations ----------------------------------------------------

#' Calculate the effect
#'
#' @param meta_data data.table. The meta information about the experiment in
#' which the contrast info can be found.
#' @param main_contrast String. Which column contains the main groups you want
#' to calculate the Hedge's G effect for. Every permutation of the groups
#' will be tested if `contrast_list` is `NULL`.
#' @param normalised_counts Numeric Matrix. The normalized count matrix.
#' @param contrast_list String vector or NULL. Optional string vector of
#' contrast formatted as `"contrast1-contrast2"`. Default NULL will create all
#' contrasts automatically.
#' @param small_sample_correction Can be NULL (automatic determination if a
#' small sample size correction should be applied) or Boolean.
#' @param parallel Boolean. Shall parallelisation be used across the contrasts.
#' @param no_cores Optional integer. Number of cores to use for parallelisation.
#' If `NULL` will default to half the detected cores.
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @returns A data.table with the effect sizes and standard errors based on the
#' Hedge's G effect size for the groups.
#'
#' @export
#'
#' @import data.table
#' @importFrom magrittr `%>%`
hedges_g_dge <- function(
  meta_data,
  main_contrast,
  normalised_counts,
  contrast_list = NULL,
  small_sample_correction = NULL,
  parallel = TRUE,
  no_cores = NULL,
  .verbose = TRUE
) {
  # checks

  checkmate::assertDataFrame(meta_data)
  checkmate::qassert(main_contrast, "S1")
  checkmate::assertClass(normalised_counts, "matrix")
  checkmate::assertNames(
    names(meta_data),
    must.include = main_contrast
  )
  checkmate::qassert(contrast_list, c("S+", "0"))
  checkmate::assert(all(grepl("-", contrast_list)))
  checkmate::qassert(parallel, "B1")
  checkmate::qassert(no_cores, c("I1", "0"))

  # function

  ## create automatic contrasts or keep the provided ones
  if (is.null(contrast_list)) {
    groups <- as.character(unique(meta_data[[main_contrast]]))
    combinations_to_test <- combn(
      x = groups,
      m = 2,
      FUN = function(x) {
        c(x[[1]], x[[2]])
      },
      simplify = FALSE
    )
  } else {
    combinations_to_test <- purrr::map_vec(
      contrast_list,
      ~ {
        x = stringr::str_split(.x, pattern = "-")
      }
    )
  }

  # parallelisation stuff
  if (!parallel) {
    mirai::daemons(1)
  } else {
    sessions <- ifelse(is.null(no_cores), get_cores(), no_cores)
    mirai::daemons(sessions)
  }

  res <- mirai::mirai_map(
    combinations_to_test,
    \(combination, meta_data, normalised_counts) {
      grpA <- meta_data[
        eval(parse(
          text = paste0(main_contrast, " == '", combination[[1]], "'")
        )),
        sample_id
      ]
      grpB <- meta_data[
        eval(parse(
          text = paste0(main_contrast, " == '", combination[[2]], "'")
        )),
        sample_id
      ]

      mat_a <- t(normalised_counts[, grpA])
      mat_b <- t(normalised_counts[, grpB])

      hedges_g_effect <- calculate_effect_size(
        mat_a = mat_a,
        mat_b = mat_b,
        small_sample_correction = small_sample_correction,
        .verbose = .verbose
      ) %>%
        data.table::setDT() %>%
        .[, `:=`(
          gene_id = colnames(mat_a),
          combination = paste(combination[[1]], combination[[2]], sep = "_vs_")
        )]

      hedges_g_effect
    },
    .args = list(meta_data = meta_data, normalised_counts = normalised_counts)
  ) %>%
    data.table::rbindlist()

  mirai::daemons(0)

  return(res)
}
