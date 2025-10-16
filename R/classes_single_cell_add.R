# additional single cell classes and methods -----------------------------------

## general generics ------------------------------------------------------------

#' Get the ready obs data from various sub method
#'
#' @description
#' Helper method that creates data.tables with cell indices which were used
#' in the given analysis + the values that are to be added to the obs table
#' in the DuckDB.
#'
#' @param x An object to set gene mapping for
#' @param ... Other parameters
#'
#' @returns Returns a data.table with a cell_idx column for the cells included
#' in the analysis and additional columns to be added to the obs table.
#'
#' @export
get_obs_data <- function(x, ...) {
  UseMethod("get_obs_data")
}

# methods ----------------------------------------------------------------------

## gene proportion analysis ----------------------------------------------------

#' @rdname get_obs_data
#'
#' @export
get_obs_data.sc_proportion_res <- function(x, ...) {
  # checks
  checkmate::assertClass(x, "sc_proportion_res")

  # function body
  obs_dt <- data.table::as.data.table(unclass(x))
  obs_dt[, cell_idx := (attr(x, "cell_indices") + 1)] # was zero indexed

  return(obs_dt)
}
