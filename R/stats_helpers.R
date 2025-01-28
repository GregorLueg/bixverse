#' Calculates a harmonic sum normalised between 0 to 1.
#'
#' @description
#' The function takes in a vector of scores between 0 and 1 and calculates a harmonic sum,
#' based on the approach OpenTargets takes to do their gene - disease evidence scores, see:
#' https://platform-docs.opentargets.org/associations.
#'
#' @param x Numeric vector. Needs to be between 0 and 1.
#'
#' @return Harmonic, normalised sum of the provided scores.
#'
#' @export
OT_harmonic_score <- function(x) {
  # Checks
  checkmate::qassert(x, "R+[0,1]")
  # Function body - finally this is fast...
  x <- sort(x, decreasing = T)
  len.x <- length(x)
  denominator <- seq_len(len.x) ^ 2
  harmonic_sum <- x / denominator
  max_harmonic_sum <- rep(1, len.x) / denominator
  sum(harmonic_sum) / sum(max_harmonic_sum)
}
