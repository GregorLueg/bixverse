# helpers ----------------------------------------------------------------------

## gene normalisations ---------------------------------------------------------

#' TPM calculation
#'
#' @param counts Numeric matrix. Count matrix (gene x sample)
#' @param gene_lengths Named vector. Named vector with gene lengths.
#'
#' @returns TPM-normalised matrix.
#'
#' @export
#'
#' @examples
#' # TPM over synthetic counts with flat 2kb gene lengths
#' syn <- synthetic_bulk_cor_matrix()
#' gene_lengths <- stats::setNames(
#'   rep(2000, nrow(syn$counts)),
#'   rownames(syn$counts)
#' )
#' tpm <- calculate_tpm(syn$counts, gene_lengths)
#' colSums(tpm)[1:3]
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
#'
#' @export
#'
#' @examples
#' # RPKM over synthetic counts with flat 2kb gene lengths
#' syn <- synthetic_bulk_cor_matrix()
#' gene_lengths <- stats::setNames(
#'   rep(2000, nrow(syn$counts)),
#'   rownames(syn$counts)
#' )
#' rpkm <- calculate_rpkm(syn$counts, gene_lengths)
#' rpkm[1:3, 1:3]
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
#' BulkCoExp object.
#' @param species String. One of `c("human", "mouse", "rat")`.
#' @param ... Additional parameters passed to methods.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # median transcript length per Ensembl gene, queried from Ensembl
#' counts <- matrix(
#'   1:4,
#'   nrow = 2,
#'   dimnames = list(c("ENSG00000141510", "ENSG00000012048"), c("s1", "s2"))
#' )
#' get_gene_lengths(counts, species = "human")
#' }
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
#'
#' @examples
#' \dontrun{
#' # the matrix method, dispatched on Ensembl rownames
#' counts <- matrix(
#'   1:4,
#'   nrow = 2,
#'   dimnames = list(c("ENSG00000141510", "ENSG00000012048"), c("s1", "s2"))
#' )
#' get_gene_lengths.matrix(counts, species = "human")
#' }
#'
#' @keywords internal
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

#' @rdname get_gene_lengths
S7::method(get_gene_lengths, BulkCoExp) <- function(
  x,
  species = c("human", "mouse", "rat"),
  ...
) {
  species <- match.arg(species)

  # checks
  checkmate::assertClass(x, "bixverse::BulkCoExp")
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
#'
#' @keywords internal
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

#' Build pairwise limma contrasts
#'
#' @description
#' Stands in for `limma::makeContrasts()` for the two cases bixverse needs:
#' every pairwise difference between the levels of the main contrast, or the
#' differences given as `"a-b"` strings.
#'
#' @param coef_names String vector. The column names of the design matrix.
#' @param contrast_grps String vector. The levels of the main contrast. Only
#' used if `contrast_list` is `NULL`.
#' @param contrast_list Optional string vector of the form `"a-b"`. If `NULL`,
#' all pairwise contrasts between `contrast_grps` are built, in design column
#' order.
#'
#' @returns A named list of numeric contrast vectors, one entry per design
#' column. Names are the contrasts with `-` replaced by `_vs_`.
#'
#' @keywords internal
build_limma_contrasts <- function(
  coef_names,
  contrast_grps,
  contrast_list = NULL
) {
  # checks
  checkmate::qassert(coef_names, "S+")
  checkmate::qassert(contrast_grps, "S+")
  checkmate::qassert(contrast_list, c("S+", "0"))

  if (is.null(contrast_list)) {
    grps <- coef_names[coef_names %in% contrast_grps]
    contrast_list <- utils::combn(grps, 2, FUN = \(x) {
      paste0(x[[1]], "-", x[[2]])
    })
  }

  pairs <- strsplit(contrast_list, "-", fixed = TRUE)
  valid <- purrr::map_lgl(pairs, \(p) {
    length(p) == 2L && all(p %in% coef_names)
  })
  if (!all(valid)) {
    stop(sprintf(
      "Contrasts need the form `a-b` with both levels in the design: %s",
      paste(contrast_list[!valid], collapse = ", ")
    ))
  }

  res <- purrr::map(pairs, \(p) {
    contrast <- stats::setNames(numeric(length(coef_names)), coef_names)
    contrast[p[[1]]] <- 1
    contrast[p[[2]]] <- -1
    contrast
  })
  names(res) <- gsub("-", "_vs_", contrast_list)

  res
}

# dge functions ----------------------------------------------------------------

## limma voom dge  -------------------------------------------------------------

#' Wrapper for a Limma Voom analysis
#'
#' @description
#' Runs the limma-voom workflow (`calcNormFactors()` -> `voomLmFit()` ->
#' `contrasts.fit()` -> `eBayes()` -> `topTable()`) in Rust via the `edge-rs`
#' crate, gated against limma 3.66.0. The design is `~ 0 + main_contrast +
#' co_variates` and every requested contrast is tested separately. limma and
#' edgeR are not needed.
#'
#' @param meta_data data.table. The meta information about the experiment in
#' which the contrast info (and potential co-variates) can be found. Rows need
#' to be in the same order as the columns of `counts`.
#' @param main_contrast String. Which column contains the main groups you want
#' to test differential gene expression with the Limma-Voom workflow for.
#' @param counts Numeric matrix. Raw counts of genes x samples, with gene
#' identifiers as row names.
#' @param contrast_list String vector or NULL. Optional string vector of
#' contrast formatted as `"contrast1-contrast2"`. Default NULL will create all
#' contrasts automatically.
#' @param co_variates String or NULL. Optional co-variates you wish to consider
#' during model fitting.
#' @param limma_params List. The limma parameters, see
#' [bixverse::params_limma_voom()].
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @returns A data.table with the columns of limma's `topTable(confint = TRUE)`
#' (`gene_id`, `logFC`, `CI.L`, `CI.R`, `AveExpr`, `t`, `P.Value`,
#' `adj.P.Val`, `B`) plus `contrast`, sorted by p-value within each contrast.
#'
#' @export
#'
#' @import data.table
#' @importFrom magrittr %>%
#'
#' @references Law, et al., Genome Biol, 2014
#'
#' @examples
#' # voom fit and topTable results for the single case vs control contrast
#' syn <- synthetic_bulk_cor_matrix()
#' meta <- data.table::data.table(
#'   sample_id = colnames(syn$counts),
#'   case_control = rep(c("case", "control"), each = 50)
#' )
#' res <- run_limma_voom(
#'   meta_data = meta,
#'   main_contrast = "case_control",
#'   counts = syn$counts,
#'   .verbose = FALSE
#' )
#' head(res)
run_limma_voom <- function(
  meta_data,
  main_contrast,
  counts,
  contrast_list = NULL,
  co_variates = NULL,
  limma_params = params_limma_voom(),
  .verbose = TRUE
) {
  variables <- c(main_contrast, co_variates)
  # checks
  checkmate::assertDataFrame(meta_data)
  checkmate::qassert(main_contrast, "S1")
  checkmate::assertMatrix(
    counts,
    mode = "numeric",
    ncols = nrow(meta_data),
    row.names = "named"
  )
  checkmate::qassert(co_variates, c("S+", "0"))
  checkmate::assertNames(
    names(meta_data),
    must.include = variables
  )
  checkmate::qassert(contrast_list, c("S+", "0"))
  assertLimmaVoomParams(limma_params)
  checkmate::qassert(.verbose, "B1")

  # copy, so the caller's table is not modified by reference
  meta_data <- data.table::copy(data.table::as.data.table(meta_data))
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

  limma_contrasts <- build_limma_contrasts(
    coef_names = colnames(model_matrix),
    contrast_grps = as.character(unique(meta_data[[main_contrast]])),
    contrast_list = contrast_list
  )

  storage.mode(counts) <- "double"

  all_dge_res <- purrr::imap(limma_contrasts, \(contrast, contrast_name) {
    tested <- .resolve_tested(design = model_matrix, contrast = contrast)
    res <- rs_limma_voom(
      counts = counts,
      design = model_matrix,
      limma_params = c(limma_params, tested)
    )

    data.table::data.table(
      gene_id = rownames(counts)[res$features_to_keep],
      logFC = res$log_fc,
      CI.L = res$ci_lower,
      CI.R = res$ci_upper,
      AveExpr = res$ave_expr,
      t = res$t_stat,
      P.Value = res$p_values,
      adj.P.Val = res$fdr,
      B = res$b_stat,
      contrast = contrast_name
    ) %>%
      data.table::setorderv("P.Value")
  }) %>%
    data.table::rbindlist()

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
#' @param normalised_counts Numeric Matrix. The normalised count matrix.
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
#' @importFrom magrittr %>%
#'
#' @examples
#' # Hedge's G on log CPM counts for the case vs control contrast
#' syn <- synthetic_bulk_cor_matrix()
#' meta <- data.table::data.table(
#'   sample_id = colnames(syn$counts),
#'   case_control = rep(c("case", "control"), each = 50)
#' )
#' norm_counts <- rs_cpm(syn$counts, lib_size = NULL, log = TRUE,
#'   prior_count = 2)
#' res <- hedges_g_dge(
#'   meta_data = meta,
#'   main_contrast = "case_control",
#'   normalised_counts = norm_counts,
#'   .verbose = FALSE
#' )
#' head(res)
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
        x <- stringr::str_split(.x, pattern = "-")
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
