# methods for meta cells and their processing/metrics --------------------------

## metrics ---------------------------------------------------------------------

### meta cell purity -----------------------------------------------------------

#' Calculate meta cell purity
#'
#' @description
#' A potential metric to see how well the meta cells are aggregated is their
#' cell type purity. This helper function helps to plot the meta-cell purity
#' based on annotated cell types. These can be also just unsupervised
#' memberships to graph-based clustering, etc.
#'
#' @param object `meta_cells` class.
#' @param original_cell_type Character vector. The original cell type
#' annotations. The indices need to match with the original cell indices used
#' to generate the meta-cells! (1-indexed)
#'
#' @returns The `meta_cells` with an added columns to the observation table
#' with the purity measures
#'
#' @export
calc_meta_cell_purity <- S7::new_generic(
  name = "calc_meta_cell_purity",
  dispatch_args = "object",
  fun = function(
    object,
    original_cell_type
  ) {
    S7::S7_dispatch()
  }
)

#' @method calc_meta_cell_purity meta_cells
S7::method(calc_meta_cell_purity, meta_cells) <- function(
  object,
  original_cell_type
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, meta_cells))
  checkmate::qassert(original_cell_type, "S+")

  # calculate purity
  purity <- purrr::map_dbl(
    object[[]]$original_cell_idx,
    function(idx) {
      types <- original_cell_type[idx]
      max(table(types)) / length(types)
    }
  )

  object[["mc_purity"]] <- purity

  return(object)
}
