# checks -----------------------------------------------------------------------

#' Check correlation graph parameters
#'
#' @description Checkmate extension for checking the graph parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @export
checkCorGraphParams <- function(x) {
  # Checkmate extension
  res <- checkmate::checkList(x)
  if (!isTRUE(res))
    return(res)
  res <- checkmate::checkNames(
    names(x),
    must.include = c(
      "kernel_bandwidth",
      "min_affinity",
      "min_cor",
      "fdr_threshold"
    )
  )
  if (!isTRUE(res))
    return(res)
  res <- purrr::map_lgl(x, .f = checkmate::qtest, rules = "R1[0, 1]")
  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(
      sprintf(
        "The following element %s in graph params does not conform the expected format. \
        Needs to be doubles between 0 and 1.",
        broken_elem
      )
    )
  }
  return(TRUE)
}

#' Check resolution graph parameters
#'
#' @description Checkmate extension for checking the resolution parameters.
#'
#' @param x The list to check/assert
#'
#' @return \code{TRUE} if the check was successful, otherwise an error message.
#'
#' @export
checkCorResParams <- function(x) {
  # Checkmate extension
  res <- checkmate::checkList(x)
  if (!isTRUE(res))
    return(res)
  res <- checkmate::checkNames(names(x), must.include = c("min_res", "max_res", "number_res"))
  if (!isTRUE(res))
    return(res)
  rules <- setNames(c("R1", "R1", "I1"), c("min_res", "max_res", "number_res"))
  res <- purrr::imap_lgl(x, \(x, name) {
    checkmate::qtest(x, rules[name])
  })
  if (!isTRUE(all(res))) {
    broken_elem <- names(res)[which(!res)][1]
    return(
      sprintf(
        "The following element %s in resolution params does not conform the expected format. \
        min_res and max_res need to be doubles and number res needs to be an integer.",
        broken_elem
      )
    )
  }
  return(TRUE)
}

# asserts ----------------------------------------------------------------------

#' Assert correlation graph parameters
#'
#' @description Checkmate extension for asserting the graph parameters.
#'
#' @inheritParams checkCorGraphParams
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @export
assertCorGraphParams <- checkmate::makeAssertionFunction(checkCorGraphParams)


#' Assert resolution graph parameters
#'
#' @description Checkmate extension for asserting the resolution parameters.
#'
#' @inheritParams checkCorResParams
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is successful.
#'
#' @export
assertCorResParams <- checkmate::makeAssertionFunction(checkCorResParams)

