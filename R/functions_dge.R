# differential expression ------------------------------------------------------

## helpers ---------------------------------------------------------------------

#' Resolve what an edgeR-style test is testing
#'
#' @description
#' Turns the user-facing `coef` / `contrast` pair into the elements the Rust
#' side reads off the parameter list. `coef` accepts 1-based column positions
#' or column names of the design and is shifted to the 0-based indices Rust
#' wants. `contrast` is flattened column-major with its column count alongside.
#'
#' Supplying neither tests the last column of the design, which is what edgeR
#' and limma default to.
#'
#' @param design Numeric matrix. The design matrix of samples x coefficients.
#' @param coef Optional integer or character. The coefficient(s) to drop from
#' the null model.
#' @param contrast Optional numeric vector or matrix. Weights over the
#' coefficients, one entry (or row) per column of `design`.
#'
#' @returns A named list holding either `coef` or `contrast` plus
#' `n_contrasts`.
#'
#' @keywords internal
.resolve_tested <- function(design, coef = NULL, contrast = NULL) {
  # checks
  checkmate::assertMatrix(design, mode = "numeric")
  checkmate::qassert(coef, c("0", "I+", "S+"))

  if (!is.null(coef) && !is.null(contrast)) {
    stop("Supply either `coef` or `contrast`, not both.")
  }

  if (!is.null(contrast)) {
    checkmate::assert(
      checkmate::checkNumeric(contrast, any.missing = FALSE),
      checkmate::checkMatrix(contrast, mode = "numeric", any.missing = FALSE)
    )
    contrast <- as.matrix(contrast)
    if (nrow(contrast) != ncol(design)) {
      stop(sprintf(
        "The contrast has %i rows against %i design columns.",
        nrow(contrast),
        ncol(design)
      ))
    }
    # as.vector() on a matrix is column-major, which is the layout Rust reads
    return(list(
      contrast = as.vector(contrast),
      n_contrasts = ncol(contrast)
    ))
  }

  if (is.null(coef)) {
    coef <- ncol(design)
  }

  if (is.character(coef)) {
    checkmate::assertSubset(coef, colnames(design))
    coef <- match(coef, colnames(design))
  }

  checkmate::assertIntegerish(
    coef,
    lower = 1L,
    upper = ncol(design),
    any.missing = FALSE
  )

  # Rust reads these 0-indexed
  list(coef = as.integer(coef) - 1L)
}

## edgeR quasi-likelihood ------------------------------------------------------

#' Run the edgeR quasi-likelihood workflow
#'
#' @description
#' Runs `filterByExpr()` -> `calcNormFactors()` -> `glmQLFit()` ->
#' `glmQLFTest()` in Rust via the `edge-rs` crate, gated against edgeR 4.8.2.
#'
#' The tested axis does not have to be genes. Anything with a counts matrix of
#' features x samples goes through here, which is why the Milo neighbourhood
#' test calls the same function with `filter = FALSE`.
#'
#' By default this skips `estimateDisp()` and lets the fit estimate its own
#' dispersion from the most abundant features. That is edgeR 4's own
#' recommendation and where most of the runtime went. Set `legacy = TRUE` in
#' [bixverse::params_edger_ql()] for the pre-4.0 pipeline, which does run
#' `estimateDisp()`.
#'
#' @param counts Numeric matrix. Raw counts of features x samples, with
#' rownames. Must not be normalised or log-transformed.
#' @param design Numeric matrix. The design matrix of samples x coefficients,
#' including the intercept. Usually the output of [stats::model.matrix()].
#' Needs at least two columns, since the null model has to retain one.
#' @param coef Optional integer or character. Which coefficient(s) of `design`
#' to drop from the null model, given as 1-based column positions or column
#' names. Defaults to the last column, as edgeR does.
#' @param contrast Optional numeric vector or matrix. Weights over the
#' coefficients, one entry (or row) per column of `design`. Mutually exclusive
#' with `coef`.
#' @param edger_params A list, see [bixverse::params_edger_ql()]. The list has
#' the following parameters:
#' \itemize{
#'   \item norm_method - String. Library size normalisation. One of
#'   `c("TMM", "TMMwsp", "RLE", "upperquartile", "none")`.
#'   \item filter - Boolean. Run `filterByExpr()` before fitting.
#'   \item min_mean - Numeric. Minimum mean count across samples.
#'   \item robust - Boolean. Robust empirical Bayes squeezing.
#'   \item legacy - Boolean. Take edgeR's pre-4.0 pipeline.
#' }
#'
#' @returns A data.table with one row per feature that survived the filters:
#' \itemize{
#'   \item feature_id - The rowname of `counts`.
#'   \item log_fc - Log2 fold change of the tested coefficient or contrast.
#'   \item log_cpm - Average log2 counts per million.
#'   \item f_stat - The quasi-likelihood F statistic.
#'   \item p_value - Raw p-value.
#'   \item fdr - Benjamini-Hochberg adjusted p-value.
#' }
#'
#' @references Chen, Lun and Smyth, F1000Research, 2016
#'
#' @export
run_edger_ql <- function(
  counts,
  design,
  coef = NULL,
  contrast = NULL,
  edger_params = params_edger_ql()
) {
  # checks
  checkmate::assertMatrix(
    counts,
    mode = "numeric",
    any.missing = FALSE,
    row.names = "named"
  )
  checkmate::assertMatrix(design, mode = "numeric", any.missing = FALSE)
  assertEdgeRQlParams(edger_params)

  if (ncol(design) < 2L) {
    stop(
      paste(
        "The design needs at least two columns:",
        "the null model has to retain one."
      )
    )
  }
  if (nrow(design) != ncol(counts)) {
    stop(sprintf(
      "The design has %i rows against %i samples in the counts.",
      nrow(design),
      ncol(counts)
    ))
  }
  if (any(counts < 0)) {
    stop("`counts` holds negative values. This wants raw counts.")
  }

  tested <- .resolve_tested(design = design, coef = coef, contrast = contrast)

  res <- rs_edger_ql(
    counts = counts,
    design = design,
    edger_params = c(edger_params, tested)
  )

  kept <- res$features_to_keep

  data.table::data.table(
    feature_id = rownames(counts)[kept],
    log_fc = res$log_fc,
    log_cpm = res$log_cpm,
    f_stat = res$f_stat,
    p_value = res$p_values,
    fdr = res$fdr
  )
}

## pseudobulk ------------------------------------------------------------------

#' Run the edgeR quasi-likelihood workflow on pseudo-bulked single cells
#'
#' @description
#' Sums the raw counts per sample and treats the result as a bulk experiment.
#' That is the whole method, and it is the one that holds its nominal false
#' discovery rate when the cells within a sample are not independent, which
#' they never are.
#'
#' This is a convenience wrapper: it pseudo-bulks with
#' [bixverse::get_pseudobulked_sc()] and hands the aggregate to
#' [bixverse::run_edger_ql()]. Reach for the two separately if you want to
#' inspect or reuse the aggregated matrix.
#'
#' @param object `SingleCells` or `SingleCellsSubset` class.
#' @param cell_list Named list of character vectors. The cell identifiers per
#' pseudo-bulk sample. The names become the sample identifiers, and the rows of
#' `design` must follow that order.
#' @param design Numeric matrix. The design matrix of samples x coefficients,
#' including the intercept. Rows aligned to `cell_list`.
#' @param coef Optional integer or character. See [bixverse::run_edger_ql()].
#' @param contrast Optional numeric vector or matrix. See
#' [bixverse::run_edger_ql()].
#' @param edger_params A list, see [bixverse::params_edger_ql()].
#' @param .verbose Boolean or integer. Controls verbosity and returns run
#' times. `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` ->
#' detailed verbosity.
#'
#' @returns A data.table, see [bixverse::run_edger_ql()]. The `feature_id`
#' column holds the gene identifiers.
#'
#' @references Squair, et al., Nat Commun, 2021; Chen, Lun and Smyth,
#' F1000Research, 2016
#'
#' @export
pseudobulk_dge_sc <- function(
  object,
  cell_list,
  design,
  coef = NULL,
  contrast = NULL,
  edger_params = params_edger_ql(),
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(
    S7::S7_inherits(object, SingleCells) ||
      S7::S7_inherits(object, SingleCellsSubset)
  )
  checkmate::assertList(cell_list, types = "character", names = "named")
  checkmate::assertMatrix(design, mode = "numeric", any.missing = FALSE)
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  if (nrow(design) != length(cell_list)) {
    stop(sprintf(
      "The design has %i rows against %i pseudo-bulk samples.",
      nrow(design),
      length(cell_list)
    ))
  }

  # raw counts only: a negative binomial cannot model the mean of normalised
  # counts over a group, zeros included
  aggregate <- get_pseudobulked_sc(
    object = object,
    cell_list = cell_list,
    return_format = "dense",
    assay = "raw",
    .verbose = .verbose
  )

  run_edger_ql(
    counts = t(aggregate),
    design = design,
    coef = coef,
    contrast = contrast,
    edger_params = edger_params
  )
}
