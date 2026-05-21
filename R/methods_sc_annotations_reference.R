# cell type reference and annotation methods -----------------------------------

## sc type ---------------------------------------------------------------------

#' @export
calc_sc_type_scores <- S7::new_generic(
  name = "calc_sc_type_scores",
  dispatch_args = "object",
  fun = function(
    object,
    cell_marker_list,
    sensitivity = FALSE,
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
  sensitivity = FALSE,
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
