# helpers ----------------------------------------------------------------------

## gene normalisations ---------------------------------------------------------

#' TPM calculation
#'
#' @param counts Numeric matrix. Count matrix (gene x sample)
#' @param gene_lengths Named vector. Named vector with gene lengths.
#'
#' @returns TPM-normalised matrix.
calculate_tpm <- function(counts, gene_lengths) {
  # checks
  checkmate::assertMatrix(
    counts,
    mode = "numeric",
    row.names = "named",
    col.names = "named"
  )
  checkmate::assertNumeric(gene_lengths, names = "named", any.missing = FALSE)
  checkmate::assertTRUE(all(names(gene_lengths) == rownames(counts)))

  # calculations
  rpk <- counts / (gene_lengths / 1000)
  scaling_factors <- colSums(rpk, na.rm = TRUE)
  tpm <- sweep(rpk, 2, scaling_factors, "/") * 1e6

  return(tpm)
}

#' RPKM calculation
#'
#' @param counts Numeric matrix. Count matrix (gene x sample)
#' @param gene_lengths Named vector. Named vector with gene lengths.
#'
#' @returns RPKM-normalised matrix.
calculate_rpkm <- function(counts, gene_lengths) {
  # checks
  checkmate::assertMatrix(
    counts,
    mode = "numeric",
    row.names = "named",
    col.names = "named"
  )
  checkmate::assertNumeric(gene_lengths, names = "named", any.missing = FALSE)
  checkmate::assertTRUE(all(names(gene_lengths) == rownames(counts)))

  rpk <- counts / (gene_lengths / 1000)
  total_reads_millions <- colSums(counts, na.rm = TRUE) / 1e6
  rpkm <- sweep(rpk, 2, total_reads_millions, "/")

  return(rpkm)
}

## gene lengths ----------------------------------------------------------------

#' Get the gene lengths
#'
#' @param x Object to extract gene lengths from. Can be a matrix or
#' bulk_coexp object.
#' @param species String. One of `c("human", "mouse", "rat")`.
#' @param ... Additional parameters passed to methods.
#'
#' @export
get_gene_lengths <- function(x, species = c("human", "mouse", "rat"), ...) {
  UseMethod("get_gene_lengths")
}

#' Get gene set lengths for a matrix
#'
#' @param x Numerical matrix. Assumes genes x samples as format and Ensembl
#' identifiers as gene ids.
#' @param species String. One of `c("human", "mouse", "rat")`.
#' @param ... Additional parameters. Not in use atm.
#'
#' @returns Named numeric representing the gene lengths.
#'
#' @export
get_gene_lengths.matrix <- function(
  x,
  species = c("human", "mouse", "rat"),
  ...
) {
  species <- match.arg(species)
  checkmate::assertMatrix(
    x,
    mode = "numeric",
    row.names = "named",
    col.names = "named"
  )
  checkmate::assertChoice(species, c("human", "mouse", "rat"))

  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop(paste(
      "Package 'biomaRt' is required to use this function.",
      "Please install it with: BiocManager::install('biomaRt')"
    ))
  }

  dataset <- c(
    "human" = "hsapiens_gene_ensembl",
    "mouse" = "mmusculus_gene_ensembl",
    "rat" = "rnorvegicus_gene_ensembl"
  )

  ensembl <- biomaRt::useMart("ensembl", dataset = dataset[species])
  gene_info <- biomaRt::getBM(
    attributes = c("ensembl_gene_id", "transcript_length"),
    filters = "ensembl_gene_id",
    values = rownames(x), # Fixed: was count_matrix
    mart = ensembl
  )

  # calculate median transcript length per gene
  gene_lengths_df <- aggregate(
    transcript_length ~ ensembl_gene_id,
    data = gene_info,
    FUN = median
  )

  gene_lengths <- gene_lengths_df$transcript_length[
    match(rownames(x), gene_lengths_df$ensembl_gene_id)
  ]
  names(gene_lengths) <- rownames(x) # Fixed: was count_matrix

  if (any(is.na(gene_lengths))) {
    warning("Some of the genes were not found. Using median imputation.")
    gene_lengths[is.na(gene_lengths)] <- median(gene_lengths, na.rm = TRUE)
  }

  return(gene_lengths)
}

#' @name get_gene_lengths.bulk_coexp
#'
#' @title Get gene set lengths for a bulk_coexp class
#'
#' @description
#' Get gene lengths for a bulk_coexp object by extracting the count matrix
#' and delegating to the matrix method.
#'
#' @param x An object of class `bulk_coexp`.
#' @param species String. One of `c("human", "mouse", "rat")`.
#' @param ... Additional parameters. Not in use atm.
#'
#' @returns Named numeric representing the gene lengths.
#'
#' @method get_gene_lengths bulk_coexp
S7::method(get_gene_lengths, bulk_coexp) <- function(
  x,
  species = c("human", "mouse", "rat"),
  ...
) {
  species <- match.arg(species)

  # checks
  checkmate::assertClass(x, "bixverse::bulk_coexp")
  checkmate::assertChoice(species, c("human", "mouse", "rat"))

  # Fixed: condition was inverted
  if (is.null(S7::prop(x, "outputs")[["raw_counts_filtered"]])) {
    stop(paste(
      "Could not find the filtered counts in the object.",
      "Did you run qc_bulk_dge()?"
    ))
  }

  counts <- S7::prop(x, "outputs")[["raw_counts_filtered"]] # Fixed: was object

  # delegate to matrix
  get_gene_lengths(counts, species = species, ...)
}

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
      "_{2,}",
      "_",
      gsub(
        "\\.",
        "_",
        make.names(gsub("[[:punct:]&&[^_]]", "", x))
      )
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
prep_limma_contrasts <- function(limma_fit, contrast_list) {
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

## limma voom dge  -------------------------------------------------------------

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
#' @param contrast_list String vector or NULL. Optional string vector of
#' contrast formatted as `"contrast1-contrast2"`. Default NULL will create all
#' contrasts automatically.
#' @param co_variates String or NULL. Optional co-variates you wish to consider
#' during model fitting.
#' @param quantile_norm Boolean. Shall the counts be also quantile-normalised.
#' Defaults to `FALSE`.
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
  contrast_list = NULL,
  co_variates = NULL,
  quantile_norm = FALSE,
  ...,
  .verbose = TRUE
) {
  variables <- c(main_contrast, co_variates)
  # Checks
  checkmate::assertDataFrame(meta_data)
  checkmate::qassert(main_contrast, "S1")
  checkmate::checkClass(dge_list, "DGEList")
  checkmate::qassert(co_variates, c("S+", "0"))
  checkmate::assertNames(
    names(meta_data),
    must.include = variables
  )
  checkmate::qassert(.verbose, "B1")
  checkmate::qassert(contrast_list, c("S+", "0"))
  checkmate::assert(all(grepl("-", contrast_list)))

  # deal with metadata columns
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

  voom_obj <- limma::voom(
    counts = dge_list,
    design = model_matrix,
    normalize.method = ifelse(quantile_norm, "quantile", "none"),
    plot = FALSE,
    ...
  )
  limma_fit <- limma::lmFit(voom_obj, model_matrix)

  contrast_grps <- unique(meta_data[[main_contrast]])

  if (is.null(contrast_list)) {
    limma_contrasts <- all_limma_contrasts(
      limma_fit = limma_fit,
      contrast_grps = contrast_grps
    )
  } else {
    limma_contrasts <- prep_limma_contrasts(
      limma_fit = limma_fit,
      contrast_list = contrast_list
    )
  }

  all_dge_res <- purrr::map(
    limma_contrasts,
    \(contrast_obj) {
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
    }
  ) %>%
    rbindlist()

  return(all_dge_res)
}

## effect size calculations ----------------------------------------------------

#' Calculate the effect size
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

  res <- purrr::map(
    combinations_to_test,
    \(combination) {
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
        small_sample_correction = NULL,
        .verbose = .verbose
      ) %>%
        data.table::setDT() %>%
        .[, `:=`(
          gene_id = colnames(mat_a),
          combination = paste(combination[[1]], combination[[2]], sep = "_vs_")
        )]

      hedges_g_effect
    }
  ) %>%
    data.table::rbindlist()

  return(res)
}
