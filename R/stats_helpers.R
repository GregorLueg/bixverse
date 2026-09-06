# various statistical functions ------------------------------------------------

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
#' @returns Harmonic, normalised sum of the provided scores.
#'
#' @export
#'
#' @examples
#' # harmonic sum over four evidence scores
#' ot_harmonic_score(c(1, 0.8, 0.5, 0.1))
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
#' @returns x, robustly scaled.
#'
#' @export
#'
#' @examples
#' # median-centred, IQR-scaled vector
#' set.seed(123)
#' head(robust_scale(rnorm(10)))
robust_scale <- function(x) {
  # Checks
  checkmate::qassert(x, "r+")
  (x - median(x, na.rm = T)) / IQR(x, na.rm = T)
}


#' Calculate the Hedge's G effect between two matrices
#'
#' @description
#' This function takes two matrices in and calculate on a per column basis the
#' Hedge's G effect size and the standard error. These results can be
#' subsequently used for meta-analyses or other approaches.
#'
#' @param mat_a Numerical matrix. Contains the values for group a. Assumes that
#' rows = samples, and columns = features.
#' @param mat_b Numerical matrix. Contains the values for group b.
#' @param small_sample_correction Can be NULL (automatic determination if a
#' small sample size correction should be applied) or Boolean.
#' @param .verbose Boolean that controls verbosity of the function.
#'
#' @returns x, robustly scaled.
#'
#' @export
#'
#' @examples
#' # Hedge's G between two groups of ten samples
#' set.seed(42)
#' mat_a <- matrix(rnorm(100), nrow = 10, ncol = 10)
#' mat_b <- matrix(rnorm(100, mean = 1), nrow = 10, ncol = 10)
#' colnames(mat_a) <- colnames(mat_b) <- sprintf("gene_%i", 1:10)
#' res <- calculate_effect_size(mat_a, mat_b, .verbose = FALSE)
#' head(res$effect_sizes)
calculate_effect_size <- function(
  mat_a,
  mat_b,
  small_sample_correction = NULL,
  .verbose = TRUE
) {
  # Checks
  checkmate::assertMatrix(mat_a, mode = "numeric", min.rows = 3L, min.cols = 1L)
  checkmate::assertMatrix(mat_b, mode = "numeric", min.rows = 3L, min.cols = 1L)
  checkmate::qassert(small_sample_correction, c("B1", "0"))
  # Function
  intersecting_features <- intersect(colnames(mat_a), colnames(mat_b))
  mat_a <- mat_a[, intersecting_features]
  mat_b <- mat_b[, intersecting_features]

  message_text <- if (!is.null(small_sample_correction)) {
    sprintf(
      "Using user-specified choice for small sample correction. Correction is set to %s",
      small_sample_correction
    )
  } else {
    total_n <- nrow(mat_a) + nrow(mat_b)
    ifelse(
      total_n <= 50,
      "Less than 50 samples identified. Applying small sample correction.",
      "More than 50 samples identified. No small sample correction applied."
    )
  }

  if (.verbose) {
    message(message_text)
  }

  small_sample_correction <- if (is.null(small_sample_correction)) {
    total_n <= 50
  } else {
    small_sample_correction
  }

  # TO DO: implement Glass effect size estimation

  results <- rs_hedges_g(
    mat_a = mat_a,
    mat_b = mat_b,
    small_sample_correction = small_sample_correction
  )

  results
}

#' Find a threshold via the Otsu method
#'
#' @param x Numerical vector. The vector for which to find the threshold via
#' the Otsu method
#' @param bins Integer. Number of bins to use.
#'
#' @returns The threshold.
#'
#' @export
#'
#' @examples
#' # threshold separating two Gaussian modes
#' set.seed(42)
#' find_threshold_otsu(c(rnorm(100), rnorm(100, mean = 5)))
find_threshold_otsu <- function(x, bins = 100L) {
  checkmate::qassert(x, "N+")
  checkmate::qassert(bins, "I+")
  checkmate::assertTRUE(length(x) > length(bins))

  rs_sc_otsu_method(scores = x, bins = bins)
}
