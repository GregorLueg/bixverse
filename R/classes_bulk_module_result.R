# BulkModuleResult -------------------------------------------------------------

## constructor -----------------------------------------------------------------

#' Constructor for a bulk co-expression module result
#'
#' @description
#' Uniform container for the terminal output of every bulk co-expression
#' method. Stores the module membership (gene to module_id assignment), the
#' method-specific factor matrices (gene loadings, sample activities,
#' eigengenes, dictionaries), the fit parameters, and method-specific
#' diagnostics.
#'
#' @param modules data.table. Module membership. Must contain `gene` and
#' `module_id` columns. Method-specific extras (`sign`, `stability`,
#' `loading`, ...) are allowed.
#' @param factors Named list of matrices. Method-agnostic keys shared across
#' methods: `gene_loadings`, `sample_activity`, `module_eigengenes`,
#' `dictionary`, `loadings`, `gene_to_eigengene_cor`. May be empty for
#' methods that produce no factor matrices (e.g. Leiden).
#' @param method String. One of `c("correlation-based",`
#' `"differential correlation-based", "ICA-based", "dgrdl-based",`
#' `"nmf-based")`.
#' @param params List. The `params_xxx()` list that produced the fit.
#' @param diagnostics Named list. Method-specific diagnostics (stability
#' data.table, resolution used, per-run losses, laplacians, ...).
#'
#' @returns An object of class `BulkModuleResult`.
#'
#' @keywords internal
new_bulk_module_result <- function(
  modules,
  factors = list(),
  method,
  params = list(),
  diagnostics = list()
) {
  checkmate::assertDataTable(modules)
  checkmate::assertNames(
    names(modules),
    must.include = c("gene", "module_id")
  )
  checkmate::assertList(factors, names = "named", null.ok = FALSE)
  checkmate::assertChoice(
    method,
    c(
      "correlation-based",
      "differential correlation-based",
      "ICA-based",
      "dgrdl-based",
      "nmf-based"
    )
  )
  checkmate::assertList(params)
  checkmate::assertList(diagnostics, names = "named", null.ok = FALSE)

  res <- list(
    modules = modules,
    factors = factors,
    method = method,
    params = params,
    diagnostics = diagnostics
  )
  class(res) <- "BulkModuleResult"
  res
}

## primitives ------------------------------------------------------------------

#' @method print BulkModuleResult
#'
#' @export
print.BulkModuleResult <- function(x, ...) {
  cat("BulkModuleResult\n")
  cat(sprintf("  Method:           %s\n", x$method))
  cat(sprintf("  No genes:         %d\n", length(unique(x$modules$gene))))
  cat(sprintf(
    "  No modules:       %d\n",
    length(unique(x$modules$module_id))
  ))
  cat(sprintf(
    "  Factor keys:      %s\n",
    if (length(x$factors) == 0L) {
      "<none>"
    } else {
      paste(names(x$factors), collapse = ", ")
    }
  ))
  cat(sprintf(
    "  Diagnostics keys: %s\n",
    if (length(x$diagnostics) == 0L) {
      "<none>"
    } else {
      paste(names(x$diagnostics), collapse = ", ")
    }
  ))
  invisible(x)
}

#' @export
#'
#' @keywords internal
format.BulkModuleResult <- function(x, ...) {
  sprintf(
    "BulkModuleResult (%s, %d genes, %d modules)",
    x$method,
    length(unique(x$modules$gene)),
    length(unique(x$modules$module_id))
  )
}

#' @export
#'
#' @keywords internal
dim.BulkModuleResult <- function(x) {
  c(
    length(unique(x$modules$gene)),
    length(unique(x$modules$module_id))
  )
}

## generics --------------------------------------------------------------------

#' Get the module membership from a BulkModuleResult
#'
#' @description
#' Returns the data.table of gene to module_id assignments. The exact
#' columns depend on the method that produced the result (CoReMo adds
#' `sign` and `stability`; NMF/ICA/DGRDL add `loading` and `sign`;
#' Leiden adds only `module_id`).
#'
#' @param object A `BulkModuleResult`.
#'
#' @returns data.table with at minimum `gene` and `module_id` columns.
#'
#' @export
get_modules <- S7::new_generic(
  name = "get_modules",
  dispatch_args = "object",
  fun = function(object) {
    S7::S7_dispatch()
  }
)

#' @method get_modules BulkModuleResult
#'
#' @export
S7::method(get_modules, S7::new_S3_class("BulkModuleResult")) <-
  function(object) {
    checkmate::assertClass(object, "BulkModuleResult")
    object[["modules"]]
  }

#' Get factor matrices from a BulkModuleResult
#'
#' @description
#' Returns one factor matrix by key, or the whole named list of factor
#' matrices if `which` is `NULL`. Keys are method-specific; see the
#' `factors` argument of `new_bulk_module_result()` for the shared vocabulary.
#'
#' @param object A `BulkModuleResult`.
#' @param which String or `NULL`. Name of the factor matrix to return. If
#' `NULL`, returns the whole list.
#'
#' @returns A matrix, the named list of matrices, or `NULL` (with warning)
#' if `which` is not among the stored factor keys.
#'
#' @export
get_factors <- S7::new_generic(
  name = "get_factors",
  dispatch_args = "object",
  fun = function(object, which = NULL) {
    S7::S7_dispatch()
  }
)

#' @method get_factors BulkModuleResult
#'
#' @export
S7::method(get_factors, S7::new_S3_class("BulkModuleResult")) <-
  function(object, which = NULL) {
    checkmate::assertClass(object, "BulkModuleResult")
    checkmate::qassert(which, c("0", "S1"))
    factors <- object[["factors"]]
    if (is.null(which)) {
      return(factors)
    }
    if (!(which %in% names(factors))) {
      warning(sprintf(
        "Factor `%s` not found in this BulkModuleResult. Available: %s. Returning NULL.",
        which,
        if (length(factors) == 0L) {
          "<none>"
        } else {
          paste(names(factors), collapse = ", ")
        }
      ))
      return(NULL)
    }
    factors[[which]]
  }

#' Get diagnostics from a BulkModuleResult
#'
#' @description
#' Returns one diagnostic by key, or the whole named list if `which` is
#' `NULL`. Contents are method-specific: CoReMo stores `stability` and
#' `cluster_quality`; Leiden stores `resolution_used` and `modularity`; ICA
#' stores `ica_meta` and `stability_scores`; DGRDL stores
#' `feature_laplacian`, `sample_laplacian`, and the grid-search table; NMF
#' stores `final_loss`, `n_iter`, `converged`, and (for stabilised runs)
#' `losses`, `best_idx`.
#'
#' @param object A `BulkModuleResult`.
#' @param which String or `NULL`. Name of the diagnostic to return. If
#' `NULL`, returns the whole list.
#'
#' @returns The requested diagnostic, the named list, or `NULL` (with
#' warning) if `which` is not among the stored diagnostic keys.
#'
#' @export
get_diagnostics <- S7::new_generic(
  name = "get_diagnostics",
  dispatch_args = "object",
  fun = function(object, which = NULL) {
    S7::S7_dispatch()
  }
)

#' @method get_diagnostics BulkModuleResult
#'
#' @export
S7::method(get_diagnostics, S7::new_S3_class("BulkModuleResult")) <-
  function(object, which = NULL) {
    checkmate::assertClass(object, "BulkModuleResult")
    checkmate::qassert(which, c("0", "S1"))
    diagnostics <- object[["diagnostics"]]
    if (is.null(which)) {
      return(diagnostics)
    }
    if (!(which %in% names(diagnostics))) {
      warning(sprintf(
        "Diagnostic `%s` not found in this BulkModuleResult. Available: %s. Returning NULL.",
        which,
        if (length(diagnostics) == 0L) {
          "<none>"
        } else {
          paste(names(diagnostics), collapse = ", ")
        }
      ))
      return(NULL)
    }
    diagnostics[[which]]
  }
