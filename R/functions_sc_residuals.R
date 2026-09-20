# single cell residuals --------------------------------------------------------

# scTransform and analytic Pearson share everything downstream of the fit, so
# the branches `find_hvg_sc()` and `calculate_pca_sc()` take live here rather
# than being written out once per class.

## the fit class ---------------------------------------------------------------

#' Construct the fitted residual model side-car
#'
#' @description
#' Thin wrapper over the list [bixverse::rs_sc_fit_residuals()] returns, with
#' the provenance R knows about and Rust does not.
#'
#' `covariate_names` is kept separately even though the covariates carry their
#' own names: the fitted coefficients are matched to covariates by position, so
#' a reordered selection would be applied to the wrong column. Rust checks the
#' order on every use, and holding the expected order here lets R say what went
#' wrong first.
#'
#' @param res List. The result of [bixverse::rs_sc_fit_residuals()].
#' @param method String. The method that was fitted.
#' @param params List. The parameters used.
#' @param group_column String or `NULL`. The grouping column, if any.
#'
#' @returns The `ScResidualFit` object.
#'
#' @keywords internal
new_sc_residual_fit <- function(res, method, params, group_column = NULL) {
  # checks
  checkmate::assertList(res)
  checkmate::assertChoice(method, c("sctransform", "analytic_pearson"))
  checkmate::assertList(params)
  checkmate::qassert(group_column, c("S1", "0"))

  res[["params"]] <- params
  res[["group_column"]] <- if (is.null(group_column)) {
    NA_character_
  } else {
    group_column
  }
  res[["covariate_names"]] <- names(res[["covariates"]])

  class(res) <- "ScResidualFit"

  return(res)
}

## primitives ------------------------------------------------------------------

#' @export
print.ScResidualFit <- function(x, ...) {
  method <- if (identical(x$method, "sctransform")) {
    "scTransform (v2)"
  } else {
    "analytic Pearson"
  }

  cat(sprintf("Fitted residual model: %s\n", method))
  cat(sprintf(
    "  %s over %s cells\n",
    ngettext(x$n_groups, "1 group", sprintf("%i groups", x$n_groups)),
    length(x$cell_indices)
  ))
  cat(sprintf("  %s genes modelled in every group\n", length(x$genes)))

  if (!is.na(x$group_column)) {
    cat(sprintf("  grouped by: %s\n", x$group_column))
  }
  if (length(x$covariate_names) > 0) {
    cat(sprintf(
      "  covariates: %s\n",
      paste(x$covariate_names, collapse = ", ")
    ))
  }

  invisible(x)
}

## guards ----------------------------------------------------------------------

#' Fetch a fitted residual model that is safe to compute with
#'
#' @description
#' The getter warns and hands back `NULL` when nothing was fitted, which is
#' what a presence probe wants. Anything about to compute residuals needs the
#' harder version, so this errors instead.
#'
#' Also compares the fitted cell set against the current one. Rust makes the
#' same check, but by the time it fires the message is about index vectors; here
#' it can name the function to re-run.
#'
#' @param object `SingleCells` or `SingleCellsSubset` class.
#' @param cell_indices Integer. The 0-based cells about to be used.
#'
#' @returns The `ScResidualFit`.
#'
#' @keywords internal
.assert_fit_usable <- function(object, cell_indices) {
  fit <- suppressWarnings(get_residual_fit(object))

  if (is.null(fit)) {
    stop(paste(
      "No fitted residual model found on the object.",
      "Run fit_residuals_sc() first."
    ))
  }

  if (!identical(as.integer(fit$cell_indices), as.integer(cell_indices))) {
    stop(sprintf(
      paste(
        "The residual model was fitted on %i cells but %i are selected now.",
        "The cell filter moved after fitting; re-run fit_residuals_sc()."
      ),
      length(fit$cell_indices),
      length(cell_indices)
    ))
  }

  return(fit)
}

#' Restrict genes to those the fitted model covers
#'
#' @description
#' A gene one group filtered out never reaches the shared axis, so an HVG set
#' picked by another method can name genes the model has no coefficients for.
#' Rust would error naming a raw store index, which is not something you can act
#' on from R.
#'
#' @param gene_indices Integer. The 0-based genes requested.
#' @param fit `ScResidualFit` object.
#' @param .verbose Boolean or Integer. Controls verbosity.
#'
#' @returns The 0-based genes the model covers, ascending.
#'
#' @keywords internal
.residual_gene_axis <- function(gene_indices, fit, .verbose = TRUE) {
  covered <- sort(intersect(as.integer(gene_indices), as.integer(fit$genes)))
  dropped <- length(gene_indices) - length(covered)

  if (length(covered) == 0) {
    stop(paste(
      "None of the requested genes are covered by the fitted residual model.",
      "The model only covers genes detected in every group."
    ))
  }

  if (dropped > 0 && .verbose) {
    message(sprintf(
      paste(
        "Dropping %i of %i genes not covered by the residual model.",
        "They were filtered out in at least one group."
      ),
      dropped,
      length(gene_indices)
    ))
  }

  return(covered)
}

## workers ---------------------------------------------------------------------

#' Observation rows for the selected cells, in the selected order
#'
#' @description
#' The group labels and the covariates are matched to cells by position, so the
#' observation rows have to arrive in the same order as `cell_indices`.
#' `get_sc_obs()` reads from DuckDB without an `ORDER BY`, so the order it
#' returns is not something to rely on. Reordering here is cheap and removes a
#' failure that would otherwise be silent: a model fitted against shuffled
#' covariates gives plausible residuals and a wrong embedding.
#'
#' @param object `SingleCells` or `SingleCellsSubset` class.
#' @param cell_indices Integer. The 0-based cells, in the order Rust gets them.
#'
#' @returns The observation table for those cells, row `i` being
#' `cell_indices[i]`.
#'
#' @keywords internal
.residual_obs <- function(object, cell_indices) {
  obs <- get_sc_obs(object, filtered = TRUE)

  if (!"cell_idx" %in% names(obs)) {
    stop("The observation table has no `cell_idx` column to align on.")
  }

  # `cell_idx` is 1-based, the Rust indices are 0-based
  order_idx <- match(as.integer(cell_indices) + 1L, as.integer(obs$cell_idx))

  if (anyNA(order_idx)) {
    stop(paste(
      "Some selected cells are missing from the observation table.",
      "The object state and the database have drifted apart."
    ))
  }

  obs[order_idx, ]
}

#' Fit a residual model over a set of cells
#'
#' @description
#' Shared body of the `fit_residuals_sc()` methods. Resolves the parameters,
#' the grouping and the covariates, then hands over to Rust.
#'
#' @inheritParams fit_residuals_sc
#' @param cell_indices Integer. The 0-based cells to fit on.
#' @param obs data.table. The observation table for those cells, in the same
#' order.
#'
#' @returns The `ScResidualFit`.
#'
#' @keywords internal
.fit_residuals <- function(
  object,
  cell_indices,
  obs,
  method,
  group_column,
  covariate_columns,
  residual_params,
  gene_batch_size,
  seed,
  .verbose
) {
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

  if (.verbose) {
    message(sprintf(
      "Fitting %s over %i cells%s.",
      method,
      length(cell_indices),
      if (is.null(group_of_cell)) {
        ""
      } else {
        sprintf(" in %i groups", length(unique(group_of_cell)))
      }
    ))
  }

  res <- rs_sc_fit_residuals(
    f_path_gene = get_rust_count_gene_f_path(object),
    f_path_cell = get_rust_count_cell_f_path(object),
    method = method,
    cell_indices = cell_indices,
    group_of_cell = group_of_cell,
    covariates = covariates,
    params = residual_params,
    gene_batch_size = gene_batch_size,
    seed = seed,
    verbose = parse_verbosity(.verbose)
  )

  new_sc_residual_fit(
    res = res,
    method = method,
    params = residual_params,
    group_column = group_column
  )
}

#' Resolve and validate the parameters for a residual fit
#'
#' @param residual_params List or `NULL`. What the caller passed.
#' @param method String. The method being fitted.
#'
#' @returns The validated parameter list.
#'
#' @keywords internal
.resolve_residual_params <- function(residual_params, method) {
  residual_params <- if (!is.null(residual_params)) {
    residual_params
  } else if (method == "sctransform") {
    params_sc_sctransform()
  } else {
    params_sc_apr()
  }

  if (method == "sctransform") {
    assertScSctransform(residual_params)
  } else {
    assertScApr(residual_params)
  }

  return(residual_params)
}

#' Pull covariate columns out of the observation table
#'
#' @description
#' The design takes numerics only. A factor would need dummy coding, which
#' changes the rank of the design matrix, so it is refused rather than guessed
#' at. The library size is never a covariate: it enters the model as a fixed
#' offset, which is what separates v2 from v1.
#'
#' @param obs data.table. The observation table for the selected cells.
#' @param covariate_columns Character or `NULL`. Columns to use.
#' @param method String. The method being fitted.
#'
#' @returns A named list of numeric vectors, empty when there are none.
#'
#' @keywords internal
.residual_covariates <- function(obs, covariate_columns, method) {
  if (is.null(covariate_columns) || length(covariate_columns) == 0) {
    return(list())
  }

  if (method != "sctransform") {
    stop(paste(
      "Covariates are only supported for method = 'sctransform'.",
      "The analytic Pearson model has no design matrix."
    ))
  }

  missing_cols <- setdiff(covariate_columns, names(obs))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Covariate columns not found in the observation table: %s.",
      paste(missing_cols, collapse = ", ")
    ))
  }

  covariates <- lapply(covariate_columns, \(col) {
    values <- unlist(obs[[col]])
    if (is.factor(values) || is.character(values)) {
      stop(sprintf(
        paste(
          "Covariate '%s' is not numeric. Dummy code it yourself and pass the",
          "resulting columns; bixverse will not guess a contrast for you."
        ),
        col
      ))
    }
    as.numeric(values)
  })
  names(covariates) <- covariate_columns

  return(covariates)
}

#' Residual-based HVG selection
#'
#' @description
#' Shared body of the residual branch in `find_hvg_sc()`. The selection itself
#' happens in Rust, so this resolves the fit, calls over and writes the result
#' back.
#'
#' @param object `SingleCells` or `SingleCellsSubset` class.
#' @param hvg_no Integer. Variable features per group.
#' @param write_var Boolean. Write the residual variance to the variable table.
#' @param gene_batch_size Integer or `NULL`. Genes held in memory per batch.
#' @param .verbose Boolean or Integer. Controls verbosity.
#'
#' @returns The object with the HVGs set.
#'
#' @keywords internal
.find_hvg_residual <- function(
  object,
  hvg_no,
  write_var,
  gene_batch_size = NULL,
  .verbose = TRUE
) {
  cell_indices <- get_cells_to_keep(object)
  fit <- .assert_fit_usable(object, cell_indices)

  res <- rs_sc_residual_variance(
    f_path_gene = get_rust_count_gene_f_path(object),
    residual_fit = fit,
    cell_indices = cell_indices,
    n_hvg = hvg_no,
    gene_batch_size = gene_batch_size,
    verbose = parse_verbosity(.verbose)
  )

  if (write_var) {
    object <- set_sc_new_var_cols(
      object = object,
      data_list = list(
        residual_variance = .scatter_gene_values(
          values = rowMeans(res$variance),
          gene_indices = res$genes,
          n_genes = dim(object)[2]
        )
      )
    )
  }

  if (.verbose && fit$n_groups > 1L) {
    message(sprintf(
      paste(
        "Selected %i variable features: the union of the top %i of each of",
        "%i groups."
      ),
      length(res$hvg),
      hvg_no,
      fit$n_groups
    ))
  }

  # rust hands back 0-based indices, set_hvg() takes 1-based
  set_hvg(object, hvg = res$hvg + 1L)
}

#' Scatter per-gene values into a full-length vector
#'
#' @description
#' The variable table takes one row per gene in the store, but a residual model
#' only covers the genes it retained. Everything else is `NA` rather than zero,
#' which would read as a real measurement of no variance.
#'
#' @param values Numeric. One value per entry in `gene_indices`.
#' @param gene_indices Integer. The 0-based genes `values` belongs to.
#' @param n_genes Integer. Total number of genes in the store.
#'
#' @returns A numeric vector of length `n_genes`.
#'
#' @keywords internal
.scatter_gene_values <- function(values, gene_indices, n_genes) {
  # checks
  checkmate::qassert(values, "N+")
  checkmate::qassert(n_genes, "X1")
  checkmate::assertTRUE(length(values) == length(gene_indices))

  out <- rep(NA_real_, n_genes)
  out[as.integer(gene_indices) + 1L] <- values

  return(out)
}

## meta cells ------------------------------------------------------------------

# meta cells hold their counts in memory and index genes 1-based, so the three
# workers above do not transfer. The Rust side is the same code behind an
# in-memory reader.

#' Residual-based HVG selection for meta cells
#'
#' @param object `MetaCells` class.
#' @param hvg_no Integer. Variable features per group.
#' @param .verbose Boolean or Integer. Controls verbosity.
#'
#' @returns The object with the HVGs set.
#'
#' @keywords internal
.find_hvg_residual_mc <- function(object, hvg_no, .verbose = TRUE) {
  fit <- .assert_fit_usable_mc(object)

  res <- rs_mc_residual_variance(
    sparse_data = mc_counts_to_list(object = object, assay = "raw"),
    residual_fit = fit,
    n_hvg = hvg_no,
    verbose = parse_verbosity(.verbose)
  )

  if (.verbose && fit$n_groups > 1L) {
    message(sprintf(
      "Selected %i variable features across %i groups.",
      length(res$hvg),
      fit$n_groups
    ))
  }

  # meta cells index genes 1-based, unlike the on-disk classes
  set_hvg(object, hvg = res$hvg + 1L)
}

#' Residual-based PCA for meta cells
#'
#' @param object `MetaCells` class.
#' @param no_pcs Integer. Number of PCs.
#' @param pca_params List. See [bixverse::params_sc_pca()].
#' @param selected_hvg Integer. The 1-based genes to use.
#' @param seed Integer. Random seed.
#' @param .verbose Boolean or Integer. Controls verbosity.
#'
#' @returns The object with the PCA attached.
#'
#' @keywords internal
.calculate_pca_residual_mc <- function(
  object,
  no_pcs,
  pca_params,
  selected_hvg,
  seed,
  .verbose
) {
  .assert_residual_pca_params(pca_params)

  fit <- .assert_fit_usable_mc(object)

  gene_indices <- .residual_gene_axis(
    gene_indices = as.integer(selected_hvg) - 1L,
    fit = fit,
    .verbose = .verbose
  )

  if (.verbose) {
    message(sprintf(
      "Using dense SVD on %s residuals for %i genes.",
      fit$method,
      length(gene_indices)
    ))
  }

  # the full matrix goes over, since the fit is indexed against every gene
  res <- rs_mc_pca_residuals(
    sparse_data = mc_counts_to_list(object = object, assay = "raw"),
    residual_fit = fit,
    no_pcs = no_pcs,
    pca_params = pca_params,
    gene_indices = gene_indices,
    seed = seed,
    verbose = parse_verbosity(.verbose)
  )

  object <- set_pca_factors(object, res$scores, from = "residuals")
  object <- set_pca_loadings(object, res$loadings)
  object <- set_pca_singular_vals(object, res$singular_values[1:no_pcs])

  return(object)
}

#' Fetch a meta cell residual fit that is safe to compute with
#'
#' @description
#' Meta cells have no cell filter to drift, their cell set is fixed at
#' construction, so this is the presence check alone.
#'
#' @param object `MetaCells` class.
#'
#' @returns The `ScResidualFit`.
#'
#' @keywords internal
.assert_fit_usable_mc <- function(object) {
  fit <- suppressWarnings(get_residual_fit(object))

  if (is.null(fit)) {
    stop(paste(
      "No fitted residual model found on the object.",
      "Run fit_residuals_sc() first."
    ))
  }

  return(fit)
}

#' Reject PCA settings the residual path cannot honour
#'
#' @description
#' Both are refused rather than overridden. They are settings the caller passed
#' explicitly, and quietly changing them would make the recorded parameters
#' disagree with what was actually run.
#'
#' @param pca_params List. See [bixverse::params_sc_pca()].
#'
#' @returns Invisibly `TRUE`; called for the error.
#'
#' @keywords internal
.assert_residual_pca_params <- function(pca_params) {
  if (pca_params$clr) {
    stop(paste(
      "Residual PCA cannot apply the PFlogPF transformation, which belongs to",
      "the normalised layer. Pass params_sc_pca(clr = FALSE)."
    ))
  }

  if (pca_params$normalise_variance) {
    stop(paste(
      "Residual PCA cannot normalise the variance: the residuals already carry",
      "the biological signal as variance, and scaling it away is the one thing",
      "the transform exists to avoid.",
      "Pass params_sc_pca(normalise_variance = FALSE)."
    ))
  }

  invisible(TRUE)
}

#' Residual-based PCA
#'
#' @description
#' Shared body of the residual branch in `calculate_pca_sc()`. Enforces the
#' settings the Rust side refuses, since those errors name internals rather
#' than the argument to change.
#'
#' @param object `SingleCells` or `SingleCellsSubset` class.
#' @param no_pcs Integer. Number of PCs.
#' @param pca_params List. See [bixverse::params_sc_pca()].
#' @param selected_hvg Integer. The 0-based genes to use.
#' @param sparse_svd Boolean. Must be `FALSE`.
#' @param seed Integer. Random seed.
#' @param .verbose Boolean or Integer. Controls verbosity.
#'
#' @returns The object with the PCA attached.
#'
#' @keywords internal
.calculate_pca_residual <- function(
  object,
  no_pcs,
  pca_params,
  selected_hvg,
  sparse_svd,
  seed,
  .verbose
) {
  if (sparse_svd) {
    stop(paste(
      "Residual PCA has no sparse solver: a residual column is dense even",
      "where the counts are not. Pass sparse_svd = FALSE."
    ))
  }

  .assert_residual_pca_params(pca_params)

  cell_indices <- get_cells_to_keep(object)
  fit <- .assert_fit_usable(object, cell_indices)

  selected_hvg <- .residual_gene_axis(
    gene_indices = selected_hvg,
    fit = fit,
    .verbose = .verbose
  )

  # the residual path densifies every selected gene, so the footprint is the
  # full cells-by-genes matrix rather than the stored non-zeros
  n_cells <- length(cell_indices)
  gb <- n_cells * length(selected_hvg) * 8 / 1024^3
  if (gb > 8) {
    warning(sprintf(
      paste(
        "Residual PCA densifies the matrix: %i cells by %i genes is about",
        "%.1f GB. Reduce the HVG count or subset the cells if that is too",
        "much for this machine."
      ),
      n_cells,
      length(selected_hvg),
      gb
    ))
  }

  if (.verbose) {
    message(sprintf(
      "Using dense SVD on %s residuals for %i genes.",
      fit$method,
      length(selected_hvg)
    ))
  }

  res <- rs_sc_pca_residuals(
    f_path_gene = get_rust_count_gene_f_path(object),
    residual_fit = fit,
    no_pcs = no_pcs,
    pca_params = pca_params,
    cell_indices = cell_indices,
    gene_indices = selected_hvg,
    seed = seed,
    return_scaled = FALSE,
    verbose = parse_verbosity(.verbose)
  )

  object <- set_pca_factors(object, res$scores, from = "residuals")
  object <- set_pca_loadings(object, res$loadings)
  object <- set_pca_singular_vals(object, res$singular_values[1:no_pcs])

  return(object)
}
