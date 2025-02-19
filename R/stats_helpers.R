#' Calculates a harmonic sum normalised between 0 to 1.
#'
#' @description
#' The function takes in a vector of scores between 0 and 1 and calculates a
#' harmonic sum, based on the approach OpenTargets takes to do their gene -
#' disease evidence scores, see:
#' https://platform-docs.opentargets.org/associations.
#'
#' @param x Numeric vector. Needs to be between 0 and 1.
#'
#' @return Harmonic, normalised sum of the provided scores.
#'
#' @export
ot_harmonic_score <- function(x) {
  # Checks
  checkmate::qassert(x, "R+[0,1]")
  # Function body - using Rust here
  rs_ot_harmonic_sum(x)
}


#' Robust scaler.
#'
#' @description
#' Robust scaling, i.e., removes the median and scales data based on the
#' interquartile range (IQR). Useful if outliers are expected. NAs will be
#' ignored.
#'
#' @param x Numeric vector.
#'
#' @return x, robustly scaled.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' set.seed(123)
#' x = rnorm(10)
#'
#' x.scaled = robust_scaling(x)
#'
#' }
robust_scale <- function(x) {
  # Checks
  checkmate::qassert(x, "r+")
  (x - median(x, na.rm = T)) / IQR(x, na.rm = T)
}
