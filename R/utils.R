# general utils ----------------------------------------------------------------

## null coalescence ------------------------------------------------------------

#' Null coalescence
#'
#' @param a R object a
#' @param b R object b
#'
#' @returns If `a` is not `NULL`, a; otherwise b.
#'
#' @keywords internal
`%||%` <- function(a, b) {
  if (!is.null(a)) {
    return(a)
  } else {
    return(b)
  }
}
