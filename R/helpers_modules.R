# helpers related to co-expression module detections ---------------------------

## internal helpers for bulk co-expression methods -----------------------------

#' Resolve the target matrix for a BulkCoExp method
#'
#' @description
#' Returns the pre-processed matrix if present, otherwise falls back to the
#' raw data with a warning. Used at the top of every bulk co-expression
#' method to remove duplicated preamble.
#'
#' @param object The class, see [bixverse::BulkCoExp()].
#'
#' @returns Numeric matrix.
#'
#' @keywords internal
.get_bulk_target_mat <- function(object) {
  checkmate::assertClass(object, "bixverse::BulkCoExp")
  processed <- S7::prop(object, "processed_data")[["processed_data"]]
  if (purrr::is_empty(processed)) {
    warning("No pre-processed data found. Defaulting to the raw data.")

    return(S7::prop(object, "raw_data"))
  }
  processed
}

#' Guard on detection_method for a BulkCoExp method
#'
#' @description
#' Reads `detection_method` from the class params and checks it against the
#' set of methods the caller is designed to handle. Returns the resolved
#' method on success, `NULL` on mismatch (with a warning). If
#' `allow_unset = TRUE`, a `NULL` `detection_method` is treated as a first
#' invocation and returns `NA_character_` silently.
#'
#' @param object The class, see [bixverse::BulkCoExp()].
#' @param allowed Character vector. Detection-method strings the caller accepts.
#' @param method_label String. Human-readable label used in the warning
#' message (e.g. `"correlation"`, `"ICA"`, `"DGRDL"`, `"NMF"`).
#' @param allow_unset Boolean. If `TRUE`, a `NULL` `detection_method` is not
#' an error; returns `NA_character_`. Used by finalisers that set
#' `detection_method` themselves on first invocation.
#'
#' @returns One of:
#' \itemize{
#'  \item the resolved `detection_method` string (proceed)
#'  \item `NA_character_` (proceed, first invocation, only when
#'   `allow_unset = TRUE`)
#'  \item `NULL` (mismatch, caller should return object unchanged)
#' }
#'
#' @keywords internal
.assert_bulk_detection_method <- function(
  object,
  allowed,
  method_label,
  allow_unset = FALSE
) {
  checkmate::assertClass(object, "bixverse::BulkCoExp")
  checkmate::qassert(allowed, "S+")
  checkmate::qassert(method_label, "S1")
  checkmate::qassert(allow_unset, "B1")
  detection_method <- S7::prop(object, "params")[["detection_method"]]
  if (is.null(detection_method)) {
    if (allow_unset) {
      return(NA_character_)
    }
    warning(sprintf(
      "This class does not seem to be set for %s module detection. Returning class as is.",
      method_label
    ))
    return(NULL)
  }
  if (!(detection_method %in% allowed)) {
    warning(sprintf(
      "This class does not seem to be set for %s module detection. Returning class as is.",
      method_label
    ))
    return(NULL)
  }
  detection_method
}

#' Derive sparse module membership from a loading matrix
#'
#' @description
#' Turns a `gene x k` loading matrix from a matrix factorisation (ICA, NMF,
#' DGRDL) into a membership table by keeping the tails of each component's
#' loading distribution.
#'
#' Two thresholding rules, selected via `membership_params`:
#' \itemize{
#'  \item `"zscore"` - robust standardisation per component, centred on the
#'  median and scaled by the MAD, keeping `abs(z) > cutoff`. No distributional
#'  assumption beyond rough symmetry.
#'  \item `"fdr"` - two-sided p-values against a Normal null fitted per
#'  component by median and MAD, Benjamini-Hochberg adjusted, keeping
#'  `padj < fdr`.
#' }
#'
#' @param loadings Numeric matrix. `gene x k`, with row and column names.
#' @param membership_params List. See
#' [bixverse::params_module_membership()].
#'
#' @return A data.table with columns `gene`, `module_id`, `loading`, `sign` and
#' the per-component score (`z` or `padj` depending on the method). One row per
#' surviving (gene, component) pair, ordered by component then by descending
#' absolute loading.
#'
#' @keywords internal
.modules_from_loadings <- function(
  loadings,
  membership_params = params_module_membership()
) {
  module_id <- abs_loading <- NULL

  checkmate::assertMatrix(loadings, mode = "numeric", row.names = "named")
  assertModuleMembershipParams(membership_params)

  # NMF loadings are non-negative by construction, so there is no lower tail to
  # find and a two-sided rule would just shift the effective cutoff.
  one_tailed <- switch(
    membership_params$tails,
    "auto" = all(loadings >= 0),
    "upper" = TRUE,
    "both" = FALSE
  )

  per_component <- purrr::map(seq_len(ncol(loadings)), \(j) {
    vals <- loadings[, j]
    centre <- stats::median(vals)
    scale_val <- stats::mad(vals)

    # A degenerate component (all loadings identical) has no tail to speak of.
    if (!is.finite(scale_val) || scale_val <= 0) {
      return(NULL)
    }

    z <- (vals - centre) / scale_val
    stat <- if (one_tailed) z else abs(z)

    if (membership_params$method == "zscore") {
      keep <- stat > membership_params$cutoff
      score <- z[keep]
      score_name <- "z"
    } else {
      pvals <- if (one_tailed) {
        stats::pnorm(z, lower.tail = FALSE)
      } else {
        2 * stats::pnorm(abs(z), lower.tail = FALSE)
      }
      padj <- stats::p.adjust(pvals, method = "BH")
      keep <- padj < membership_params$fdr
      score <- padj[keep]
      score_name <- "padj"
    }

    if (!any(keep)) {
      return(NULL)
    }

    res <- data.table::data.table(
      gene = rownames(loadings)[keep],
      module_id = colnames(loadings)[j],
      loading = vals[keep],
      sign = ifelse(vals[keep] >= 0, "pos", "neg")
    )
    res[, (score_name) := score]
    res
  })

  res <- data.table::rbindlist(per_component)

  if (nrow(res) == 0L) {
    warning(paste(
      "No genes passed the module membership threshold.",
      "Consider loosening `membership_params`."
    ))
    return(res)
  }

  res[, abs_loading := abs(loading)]
  data.table::setorderv(res, c("module_id", "abs_loading"), c(1L, -1L))
  res[, abs_loading := NULL]

  res[]
}
