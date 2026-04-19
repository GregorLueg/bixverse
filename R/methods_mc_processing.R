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
#' @param object `MetaCells` class.
#' @param original_cell_type Character vector. The original cell type
#' annotations. The indices need to match with the original cell indices used
#' to generate the meta-cells! (1-indexed)
#'
#' @returns The `MetaCells` with an added columns to the observation table
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

#' @method calc_meta_cell_purity MetaCells
S7::method(calc_meta_cell_purity, MetaCells) <- function(
  object,
  original_cell_type
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, MetaCells))
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

## processing ------------------------------------------------------------------

### hvg ------------------------------------------------------------------------

# generic found in R/base_generics_sc.R

#' @method find_hvg_sc MetaCells
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(find_hvg_sc, MetaCells) <- function(
  object,
  hvg_no = 2000L,
  hvg_params = params_sc_hvg(),
  streaming = FALSE,
  .verbose = TRUE
) {
  checkmate::assertClass(object, "bixverse::SingleCells")
  checkmate::qassert(hvg_no, "I1")
  assertScHvg(hvg_params)
  checkmate::qassert(streaming, "B1")
  checkmate::qassert(.verbose, "B1")

  assay <- if (hvg_params$method == "vst") {
    "raw"
  } else {
    "norm"
  }

  count_list <- mc_counts_to_list(object = object, assay = assay)

  res <- with(
    hvg_params,
    rs_mc_hvg(
      sparse_data = count_list,
      hvg_method = method,
      loess_span = loess_span,
      binning = bin_method,
      n_bins = num_bin,
      clip_max = NULL,
    )
  )

  object <- set_sc_new_var_cols(object = object, data_list = res)

  hvg <- switch(
    hvg_params$method,
    "vst" = order(res$var_std, decreasing = TRUE)[1:hvg_no],
    "dispersion" = order(res$dispersion, decreasing = TRUE)[1:hvg_no],
    "meanvarbin" = order(res$dispersion_scaled, decreasing = TRUE)[1:hvg_no],
    stop("Unknown HVG method: ", hvg_params$method)
  )

  object <- set_hvg(object, hvg = hvg)

  return(object)
}
