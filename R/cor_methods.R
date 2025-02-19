# methods - simple correlations ----

#' Prepare correlation-based module detection
#'
#' @description
#' This is the generic function for doing the necessary preprocessing for
#' running independent component analysis.
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
#' of gene modules. Relevant data will be added to the `bulk_coexp` class.
#'
#' @usage ...
#'
#' @param bulk_coexp `bulk_coexp` class, see [bixverse::bulk_coexp()].
#' @param non_parametric_cors Boolean. Shall Spearman be used over Pearson
#' correlations. Defaults to `TRUE`.
#' @param .verbose Boolean. Controls verbosity of the function.
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
S7::method(cor_module_processing, bulk_coexp) <- function(bulk_coexp,
                                                          non_parametric_cors = TRUE,
                                                          .verbose = TRUE) {
  # Checks
  checkmate::assertClass(bulk_coexp, "bixverse::bulk_coexp")
  checkmate::qassert(non_parametric_cors, "B1")
  checkmate::qassert(.verbose, "B1")

  # Function body
  if (purrr::is_empty(S7::prop(bulk_coexp, "processed_data")[['processed_data']])) {
    warning("No pre-processed data found. Defaulting to the raw data")
    target_mat <- S7::prop(bulk_coexp, "raw_data")
  } else {
    target_mat <- S7::prop(bulk_coexp, "processed_data")[['processed_data']]
  }

  # Calculate the correlation matrix
  cor_matrx <- rs_cor(x = target_mat, spearman = non_parametric_cors)

  S7::prop(bulk_coexp, "processed_data")[["cor_mat"]] <- cor_matrx

  return(bulk_coexp)
}
