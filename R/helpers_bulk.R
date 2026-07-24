# internal helpers for bulk co-expression methods ------------------------------

#' Resolve the target matrix for a BulkCoExp method
#'
#' @description
#' Returns the pre-processed matrix if present, otherwise falls back to the
#' raw data with a warning. Used at the top of every bulk co-expression
#' method to remove duplicated preamble.
#'
#' @param object The class, see [bixverse::BulkCoExp()].
#' @param .verbose Boolean. If `FALSE`, suppresses the fallback warning.
#'
#' @returns Numeric matrix.
#'
#' @keywords internal
.get_bulk_target_mat <- function(object, .verbose = TRUE) {
  checkmate::assertClass(object, "bixverse::BulkCoExp")
  checkmate::qassert(.verbose, "B1")
  processed <- S7::prop(object, "processed_data")[["processed_data"]]
  if (purrr::is_empty(processed)) {
    if (.verbose) {
      warning("No pre-processed data found. Defaulting to the raw data.")
    }
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
#' Caller pattern:
#' \preformatted{
#' detection_method <- .assert_bulk_detection_method(object, "correlation-based", "correlation")
#' if (is.null(detection_method)) return(object)
#' }
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
