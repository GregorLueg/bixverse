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
#'
#' @noRd
`%||%` <- function(a, b) {
  if (!is.null(a)) {
    return(a)
  } else {
    return(b)
  }
}

## count formatting ------------------------------------------------------------

#' Format a count for a message
#'
#' @param x Integer-ish vector of counts.
#'
#' @returns The counts as a character vector, digits grouped in threes with an
#'   underscore, matching what the Rust side prints.
#'
#' @keywords internal
#'
#' @noRd
.fmt_n <- function(x) {
  format(x, big.mark = "_", trim = TRUE, scientific = FALSE)
}
