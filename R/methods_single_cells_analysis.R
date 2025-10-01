# single cell analysis methods -------------------------------------------------

## dges ------------------------------------------------------------------------

### mann whitney ---------------------------------------------------------------

find_markers_sc <- S7::new_generic(
  name = "find_markers_sc",
  dispatch_args = "object",
  fun = function(
    object,
    cells_1,
    cells_2,
    min_prop = 0.05,
    alternative = c("twosided", "greater", "less"),
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method find_markers_sc single_cell_exp
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(find_markers_sc, single_cell_exp) <- function(
  object,
  cells_1,
  cells_2,
  min_prop = 0.05,
  alternative = c("twosided", "greater", "less"),
  .verbose = TRUE
) {
  alternative <- match.arg(alternative)

  # checks
  checkmate::assertClass(object, "bixverse::single_cell_exp")
  checkmate::qassert(cells_1, "S+")
  checkmate::qassert(cells_2, "S+")
  checkmate::qassert(min_prop, "N1[0, 1]")
  checkmate::assertChoice(alternative, c("twosided", "greater", "less"))
  checkmate::qassert(.verbose, "B1")
}
