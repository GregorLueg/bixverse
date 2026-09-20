# single cell residual methods -------------------------------------------------

## fitting ---------------------------------------------------------------------

# generic found in R/base_generics_sc.R
# shared across SingleCells and SingleCellsSubset: both are backed by the same
# on-disk count files and expose the same accessors

#' @method fit_residuals_sc ScOrScSubset
S7::method(fit_residuals_sc, ScOrScSubset) <- function(
  object,
  method = c("sctransform", "analytic_pearson"),
  group_column = NULL,
  covariate_columns = NULL,
  residual_params = NULL,
  gene_batch_size = NULL,
  seed = 42L,
  .verbose = TRUE
) {
  method <- match.arg(method)

  # checks
  checkmate::assertTRUE(
    S7::S7_inherits(object, SingleCells) ||
      S7::S7_inherits(object, SingleCellsSubset)
  )
  checkmate::assertChoice(method, c("sctransform", "analytic_pearson"))
  checkmate::qassert(group_column, c("S1", "0"))
  checkmate::qassert(covariate_columns, c("S+", "0"))
  checkmate::assertList(residual_params, null.ok = TRUE)
  checkmate::qassert(gene_batch_size, c("X1[1,)", "0"))
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  cell_indices <- get_cells_to_keep(object)

  if (length(cell_indices) == 0) {
    warning(paste(
      "No cells kept in the object.",
      "Run set_cells_to_keep() first. Returning object as is."
    ))
    return(object)
  }

  obs <- .residual_obs(object, cell_indices)

  fit <- .fit_residuals(
    object = object,
    cell_indices = cell_indices,
    obs = obs,
    method = method,
    group_column = group_column,
    covariate_columns = covariate_columns,
    residual_params = residual_params,
    gene_batch_size = gene_batch_size,
    seed = seed,
    .verbose = .verbose
  )

  if (.verbose) {
    message(sprintf(
      "Fitted %i model(s) over %i genes.",
      fit$n_groups,
      length(fit$genes)
    ))
  }

  set_residual_fit(object, residual_fit = fit)
}

## meta cells ------------------------------------------------------------------

# generic found in R/base_generics_sc.R

# `analytic_pearson` is the better default here, but S7 wants the method
# formals to match the generic exactly, so the recommendation lives in the
# roxygen rather than in a different default

#' @method fit_residuals_sc MetaCells
S7::method(fit_residuals_sc, MetaCells) <- function(
  object,
  method = c("sctransform", "analytic_pearson"),
  group_column = NULL,
  covariate_columns = NULL,
  residual_params = NULL,
  gene_batch_size = NULL,
  seed = 42L,
  .verbose = TRUE
) {
  method <- match.arg(method)

  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, MetaCells))
  checkmate::assertChoice(method, c("sctransform", "analytic_pearson"))
  checkmate::qassert(group_column, c("S1", "0"))
  checkmate::qassert(covariate_columns, c("S+", "0"))
  checkmate::assertList(residual_params, null.ok = TRUE)
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  # the counts are already in memory here, so there is nothing to batch off
  # disk. Say so rather than accepting the argument and ignoring it
  if (!is.null(gene_batch_size)) {
    warning(paste(
      "`gene_batch_size` does not apply to meta cells, whose counts are",
      "already in memory. Ignoring it."
    ))
  }

  obs <- S7::prop(object, "obs_table")

  residual_params <- .resolve_residual_params(residual_params, method)

  group_of_cell <- if (is.null(group_column)) {
    NULL
  } else {
    if (!group_column %in% names(obs)) {
      stop(sprintf(
        "Column '%s' not found in the observation table.",
        group_column
      ))
    }
    as.integer(factor(unlist(obs[[group_column]]))) - 1L
  }

  covariates <- .residual_covariates(obs, covariate_columns, method)

  res <- rs_mc_fit_residuals(
    sparse_data = mc_counts_to_list(object = object, assay = "raw"),
    method = method,
    group_of_cell = group_of_cell,
    covariates = covariates,
    params = residual_params,
    seed = seed,
    verbose = parse_verbosity(.verbose)
  )

  fit <- new_sc_residual_fit(
    res = res,
    method = method,
    params = residual_params,
    group_column = group_column
  )

  if (.verbose) {
    message(sprintf(
      "Fitted %i model(s) over %i genes.",
      fit$n_groups,
      length(fit$genes)
    ))
  }

  set_residual_fit(object, residual_fit = fit)
}

## corrected counts ------------------------------------------------------------

# generic found in R/base_generics_sc.R

#' @method sct_corrected_counts_sc ScOrScSubset
S7::method(sct_corrected_counts_sc, ScOrScSubset) <- function(
  object,
  dir_out = NULL,
  build_cell_store = TRUE,
  overwrite = FALSE,
  gene_batch_size = NULL,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(
    S7::S7_inherits(object, SingleCells) ||
      S7::S7_inherits(object, SingleCellsSubset)
  )
  checkmate::qassert(dir_out, c("S1", "0"))
  checkmate::qassert(build_cell_store, "B1")
  checkmate::qassert(overwrite, "B1")
  checkmate::qassert(gene_batch_size, c("X1[1,)", "0"))
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  cell_indices <- get_cells_to_keep(object)
  fit <- .assert_fit_usable(object, cell_indices)

  if (!identical(fit$method, "sctransform")) {
    stop(paste(
      "Corrected counts need a scTransform fit; the analytic Pearson model has",
      "no corrected-count equivalent. Refit with method = 'sctransform'."
    ))
  }

  # never the data directory itself, that holds the source store
  dir_out <- if (is.null(dir_out)) {
    file.path(S7::prop(object, "dir_data"), "sct_corrected")
  } else {
    path.expand(dir_out)
  }

  f_path_gene <- file.path(dir_out, "counts_genes.bin")

  if (file.exists(f_path_gene) && !overwrite) {
    stop(sprintf(
      "A store already exists at '%s'. Pass overwrite = TRUE to replace it.",
      dir_out
    ))
  }

  dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

  res <- rs_sct_corrected_counts(
    f_path_gene = get_rust_count_gene_f_path(object),
    residual_fit = fit,
    cell_indices = cell_indices,
    f_path_out = f_path_gene,
    gene_batch_size = gene_batch_size,
    verbose = parse_verbosity(.verbose)
  )

  if (!build_cell_store) {
    if (.verbose) {
      message(sprintf(
        "Wrote %i genes by %i cells to %s.",
        res$n_genes,
        res$n_cells,
        f_path_gene
      ))
    }
    return(invisible(res))
  }

  .build_corrected_object(
    object = object,
    dir_out = dir_out,
    res = res,
    cell_indices = cell_indices,
    gene_batch_size = gene_batch_size,
    .verbose = .verbose
  )
}

#' Turn a corrected gene-major store into a `SingleCells`
#'
#' @description
#' The corrected counts come out gene-major only, and a `SingleCells` needs the
#' cell-major twin plus the database as well. The gene axis is the model's, so
#' the variable table is rebuilt from the genes that survived rather than copied
#' across: index `j` of the new store is a different gene from index `j` of the
#' old one.
#'
#' @param object `SingleCells` or `SingleCellsSubset` class.
#' @param dir_out String. Directory the store lives in.
#' @param res List. The result of [bixverse::rs_sct_corrected_counts()].
#' @param cell_indices Integer. The 0-based cells that were written.
#' @param gene_batch_size Integer or `NULL`. Genes read per batch.
#' @param .verbose Boolean or Integer. Controls verbosity.
#'
#' @returns The new `SingleCells` over the corrected counts.
#'
#' @keywords internal
.build_corrected_object <- function(
  object,
  dir_out,
  res,
  cell_indices,
  gene_batch_size,
  .verbose
) {
  if (.verbose) {
    message("Building the cell-major companion store.")
  }

  rs_sc_gene_store_to_cell_store(
    f_path_in = res$f_path,
    f_path_out = file.path(dir_out, "counts_cells.bin"),
    cells_per_phase = .CORRECTED_CELLS_PER_PHASE,
    gene_batch_size = gene_batch_size %||% .CORRECTED_GENE_BATCH,
    verbose = parse_verbosity(.verbose)
  )

  obs <- data.table::copy(.residual_obs(object, cell_indices))
  var_table <- data.table::copy(
    get_sc_var(object)[as.integer(res$genes) + 1L, ]
  )

  # both axes are renumbered by the populate helpers, and the old indices would
  # point into the source store, so they go rather than travel along wrong
  obs[, cell_idx := NULL]
  var_table[, gene_idx := NULL]

  # the populate helpers take the first column as the id, so put it there
  data.table::setcolorder(obs, c("cell_id", setdiff(names(obs), "cell_id")))
  data.table::setcolorder(
    var_table,
    c("gene_id", setdiff(names(var_table), "gene_id"))
  )

  corrected <- SingleCells(dir_data = dir_out)
  db <- get_sc_duckdb(corrected)
  db$populate_obs_from_data.table(obs_dt = obs)
  db$populate_var_from_data.table(var_dt = var_table)
  db$set_to_keep_column()

  corrected <- load_existing(corrected, .verbose = FALSE)

  if (.verbose) {
    message(sprintf(
      "Corrected store: %i cells by %i genes in %s.",
      res$n_cells,
      res$n_genes,
      dir_out
    ))
  }

  return(corrected)
}

#' Cells held in memory per phase when transposing a corrected store
#'
#' @description
#' Peak memory is roughly `cells_per_phase * mean_genes_per_cell * 12` bytes,
#' so 50k cells is a few hundred megabytes at a typical density. The source is
#' re-read once per phase, which is the trade.
#'
#' @keywords internal
.CORRECTED_CELLS_PER_PHASE <- 50000L

#' Genes read per batch when transposing a corrected store
#'
#' @keywords internal
.CORRECTED_GENE_BATCH <- 1000L
