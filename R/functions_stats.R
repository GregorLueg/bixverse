# helpers ----------------------------------------------------------------------

## topological overlap ---------------------------------------------------------

#' Calculate the TOM from a correlation matrix
#'
#' @param cor_mat Numerical matrix. The symmetric correlation matrix.
#' @param signed Boolean. Do you want to calculate the signed version. If set
#' to `FALSE`, the absolute correlation coefficients will be used.
#' @param version String. One of `c("v1", "v2")`. Defaults to `"v1"`.
#'
#' @details Calculates the topological overlap matrix from a correlation matrix.
#' The TOM is defined as:
#'
#' **Unsigned, v1:**
#'
#' \deqn{TOM_{ij} = \frac{a_{ij} + \sum_{k \neq i,j} a_{ik} a_{kj}}{\min(k_i, k_j) + 1 - a_{ij}}}
#'
#' **Signed, v1:**
#'
#' \deqn{TOM_{ij} = \frac{a_{ij} + \sum_{k \neq i,j} a_{ik} a_{kj}}{\min(k_i, k_j) + 1 - \left|a_{ij}\right|}}
#'
#' **Unsigned, v2:**
#'
#' \deqn{TOM_{ij} = 0.5 \left( a_{ij} + \frac{\sum_{k \neq i,j} a_{ik} a_{kj}}{\min(k_i, k_j) + a_{ij}} \right)}
#'
#' **Signed, v2:**
#'
#' \deqn{TOM_{ij} = 0.5 \left( a_{ij} + \frac{\sum_{k \neq i,j} a_{ik} a_{kj}}{\min(k_i, k_j) + \left|a_{ij}\right|} \right)}
#'
#' where \eqn{a_{ij}} is the affinity between nodes \eqn{i} and \eqn{j},
#' and \eqn{k_i = \sum_j a_{ij}} is the connectivity of node \eqn{i}.
#' For signed networks, connectivity is calculated as \eqn{k_i = \sum_j \left|a_{ij}\right|}.
#'
#' Version 2 uses a different normalization approach that scales the shared
#' neighbor contribution separately before combining it with the direct
#' connection strength.
#'
#' @returns A symmetric matrix of the same dimensions as `cor_mat` containing
#' the topological overlap measures.
#'
#' @export
calculate_tom <- function(cor_mat, signed, version = c("v1", "v2")) {
  version <- match.arg(version)

  # checks
  checkmate::assertMatrix(cor_mat, nrows = ncol(cor_mat), ncols = nrow(cor_mat))
  checkmate::qassert(signed, "B1")
  checkmate::assertChoice(version, c("v1", "v2"))

  # body
  if (!signed) {
    cor_mat <- abs(cor_mat)
  }

  tom_mat <- rs_tom(x = cor_mat, tom_type = version, signed = signed)

  return(tom_mat)
}


#' Calculate the TOM from an expression matrix
#'
#' @param x Numerical matrix. The expression matrix. Assumes that columns are
#' the genes, and rows the samples.
#' @param signed Boolean. Do you want to calculate the signed version. If set
#' to `FALSE`, the absolute correlation coefficients will be used.
#' @param version String. One of `c("v1", "v2")`. Defaults to `"v1"`
#' @param cor_method String. One of `c("pearson", spearman)`. Defaults to
#' `"pearson"`.
#'
#' @details Calculates the topological overlap matrix from an expression matrix.
#' It will first calculate the specified correlation matrix and then generate
#' the TOM. The TOM is defined as:
#'
#' **Unsigned, v1:**
#'
#' \deqn{TOM_{ij} = \frac{a_{ij} + \sum_{k \neq i,j} a_{ik} a_{kj}}{\min(k_i, k_j) + 1 - a_{ij}}}
#'
#' **Signed, v1:**
#'
#' \deqn{TOM_{ij} = \frac{a_{ij} + \sum_{k \neq i,j} a_{ik} a_{kj}}{\min(k_i, k_j) + 1 - \left|a_{ij}\right|}}
#'
#' **Unsigned, v2:**
#'
#' \deqn{TOM_{ij} = 0.5 \left( a_{ij} + \frac{\sum_{k \neq i,j} a_{ik} a_{kj}}{\min(k_i, k_j) + a_{ij}} \right)}
#'
#' **Signed, v2:**
#'
#' \deqn{TOM_{ij} = 0.5 \left( a_{ij} + \frac{\sum_{k \neq i,j} a_{ik} a_{kj}}{\min(k_i, k_j) + \left|a_{ij}\right|} \right)}
#'
#' where \eqn{a_{ij}} is the affinity between nodes \eqn{i} and \eqn{j},
#' and \eqn{k_i = \sum_j a_{ij}} is the connectivity of node \eqn{i}.
#' For signed networks, connectivity is calculated as \eqn{k_i = \sum_j \left|a_{ij}\right|}.
#'
#' Version 2 uses a different normalization approach that scales the shared
#' neighbor contribution separately before combining it with the direct
#' connection strength.
#'
#' @returns The topological overlap matrix.
#'
#' @export
calculate_tom_from_exp <- function(x, signed, version, cor_method) {
  # checks
  checkmate::assertMatrix(x)
  checkmate::qassert(signed, "B1")
  checkmate::assertChoice(version, c("v1", "v2"))

  spearman <- cor_method == "spearman"

  cor_mat <- rs_cor(x = x, spearman = spearman)

  # body
  if (!signed) {
    cor_mat <- abs(cor_mat)
  }

  tom_mat <- rs_tom(x = cor_mat, tom_type = version, signed = signed)

  return(tom_mat)
}

## scaling ---------------------------------------------------------------------

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
#' x <- rnorm(10)
#'
#' x.scaled <- robust_scaling(x)
#' }
robust_scale <- function(x) {
  # Checks
  checkmate::qassert(x, "r+")
  (x - median(x, na.rm = T)) / IQR(x, na.rm = T)
}

## effect sizes ----------------------------------------------------------------

#' Calculate the Hedge G effect between two matrices
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
#' @return x, robustly scaled.
#'
#' @export
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
      "Using user-specified choice for small sample correction. Correction is set to %b",
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
