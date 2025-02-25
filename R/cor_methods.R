# methods - simple correlations ----

#' Prepare correlation-based module detection
#'
#' @param bulk_coexp `bulk_coexp` class, see [bixverse::bulk_coexp()].
#' @param non_parametric_cors Boolean. Shall Spearman be used over Pearson
#' correlations. Defaults to `TRUE`.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @description
#' This function will
#'
#' @export
cor_module_processing <- S7::new_generic(
  "cor_module_processing",
  "bulk_coexp"
)

#' @name cor_module_processing
#'
#' @description
#' This function will generate necessarily data for correlation-based detection
#' of gene modules. In this case, this will be based on simple correlation
#' between the genes and the relevant data will be added to the `bulk_coexp`
#' class.
#'
#' @usage ...
#'
#' @return `bulk_coexp` with the needed data for subsequent identification of
#' correlation-based co-expression modules.
#'
#' @export
#'
#' @importFrom magrittr `%>%`
#' @importFrom magrittr `%$%`
#' @importFrom zeallot `%<-%`
#' @import data.table
#'
#' @method cor_module_processing bulk_coexp
S7::method(cor_module_processing, bulk_coexp) <- function(
    bulk_coexp,
    correlation_method = c("pearson", "spearman"),
    .verbose = TRUE)
  {
  # Checks
  checkmate::assertClass(bulk_coexp, "bixverse::bulk_coexp")
  checkmate::assertChoice(correlation_method, c("pearson", "spearman"))
  checkmate::qassert(.verbose, "B1")

  # Function body
  if (purrr::is_empty(S7::prop(bulk_coexp, "processed_data")[['processed_data']])) {
    warning("No pre-processed data found. Defaulting to the raw data")
    target_mat <- S7::prop(bulk_coexp, "raw_data")
  } else {
    target_mat <- S7::prop(bulk_coexp, "processed_data")[['processed_data']]
  }

  spearman <- if (correlation_method == 'pearson') {
    if (.verbose)
      message("Using Pearson correlations.")
    FALSE
  } else {
    if (.verbose)
      message("Using Spearman correlations.")
    TRUE
  }

  # Calculate the upper triangle of correlation matrix
  cor_diagonal <- rs_cor_upper_triangle(target_mat, spearman = spearman, shift = 0L)

  feature_names <- colnames(target_mat)
  total_len <- length(feature_names)

  feature_a <- purrr::map(1:total_len, \(idx) {
    rep(feature_names[[idx]], total_len - idx)
  })
  feature_a <- do.call(c, feature_a)

  feature_b <- purrr::map(1:total_len, \(idx) {
    remaining <- total_len - idx
    if (remaining > 0) feature_names[(idx + 1:remaining)] else character(0)
  })
  feature_b <- do.call(c, feature_b)

  correlation_res <- data.table(
    feature_a = feature_a,
    feature_b = feature_b,
    r = cor_diagonal,
    r_abs = abs(cor_diagonal)
  )

  correlation_params <- list(
    spearman = spearman,
    type = 'simple'
  )

  S7::prop(bulk_coexp, "processed_data")[["correlation_res"]] <- correlation_res
  S7::prop(bulk_coexp, "params")[["correlation_params"]] <- correlation_params

  return(bulk_coexp)
}
