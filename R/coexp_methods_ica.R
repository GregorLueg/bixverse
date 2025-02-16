# preprocessing ----

#' Prepare class for ICA
#'
#' @description
#' This is the generic function for doing the necessary preprocessing for
#' running independent component analysis.
#'
#' @export
ica_processing <- S7::new_generic(
  "ica_processing",
  "bulk_coexp"
)

#' @name ica_processing
#'
#' @description
#' This function will prepare the `bulk_coexp` for subsequent usage of the
#' contrastive ICA functions.
#'
#' @usage ...
#'
#' @param bulk_coexp `bulk_coexp` class, see [bixverse::bulk_coexp()].
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @return `bulk_coexp` with the needed data for ICA in the
#' properties of the class.
#'
#' @export
#'
#'
#' @importFrom magrittr `%>%`
#' @importFrom magrittr `%$%`
#' @import data.table
#'
#' @method ica_processing bulk_coexp
S7::method(ica_processing, bulk_coexp) <- function(bulk_coexp,
                                                   .verbose = TRUE) {
  # Checks
  checkmate::assertClass(bulk_coexp, "bixverse::bulk_coexp")
  checkmate::qassert(.verbose, "B1")

  # Function body
  if(purrr::is_empty(S7::prop(bulk_coexp_class, "processed_data")[['processed_data']])) {
    warning("No pre-processed data found. Defaulting to the raw data")
    target_mat <- S7::prop(bulk_coexp, "raw_data")
  } else {
    target_mat <- S7::prop(bulk_coexp_class, "processed_data")[['processed_data']]
  }

  # Transpose the matrix for whitening
  c(X_white, K) %<-% rs_whiten_matrix(t(target_mat))

  S7::prop(bulk_coexp, "processed_data")[["whiten_data"]] <- X_white
  S7::prop(bulk_coexp, "processed_data")[["K"]] <- K

  return(bulk_coexp)
}
