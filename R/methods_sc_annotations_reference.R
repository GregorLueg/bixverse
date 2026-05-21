# cell type reference and annotation methods -----------------------------------

## sc type ---------------------------------------------------------------------

#' Calculate ScType scores per cell
#'
#' @description
#' Implements the approach from
#'
#' @param object `SingleCells` or `SingleCellsMultiModal`
#' @param cell_marker A list, see [prepare_cell_markers()].
#' @param sensitivity Boolean. Shall shared marker genes be downweighted (like
#' in the original reference). Defaults to `TRUE`.
#' @param weight_floor Optional numeric. A value between 0 to 1 and sets the
#' floor for the weights if `sensitivity = TRUE`. If not provided, defaults to
#' `0` as in the original implementation.
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
#'
#' @returns An `ScTypeResults` results class
#'
#' @export
calc_sc_type_scores <- S7::new_generic(
  name = "calc_sc_type_scores",
  dispatch_args = "object",
  fun = function(
    object,
    cell_marker_list,
    sensitivity = TRUE,
    weight_floor = NULL,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method calc_sc_type_scores SingleCells
S7::method(calc_sc_type_scores, SingleCells) <- function(
  object,
  cell_marker_list,
  sensitivity = TRUE,
  weight_floor = NULL,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  assertCellMarkerList(cell_marker_list)
  checkmate::qassert(.verbose, c("I1", "B1"))

  res <- rs_sc_type(
    f_path = get_rust_count_gene_f_path(object),
    cell_indices = get_cells_to_keep(object),
    cell_markers = cell_marker_list,
    sensitivity = sensitivity,
    weight_floor = weight_floor,
    verbose = parse_verbosity(.verbose)
  )

  class(res) <- "ScTypeResults"

  return(res)
}
