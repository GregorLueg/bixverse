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
  checkmate::assertTRUE(S7::S7_inherits(object, MetaCells))
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
      clip_max = NULL
    )
  )

  object@var_table[, names(res) := res]

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

### pca ------------------------------------------------------------------------

# generic found in R/base_generics_sc.R

#' @method calculate_pca_sc MetaCells
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(calculate_pca_sc, MetaCells) <- function(
  object,
  no_pcs,
  randomised_svd = TRUE,
  sparse_svd = FALSE,
  hvg = NULL,
  seed = 42L,
  .verbose = TRUE
) {
  checkmate::assertTRUE(S7::S7_inherits(object, MetaCells))
  checkmate::qassert(no_pcs, "I1")
  checkmate::qassert(randomised_svd, "B1")
  checkmate::qassert(sparse_svd, "B1")
  checkmate::qassert(hvg, c("I+", "0"))
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, "B1")

  if ((length(get_hvg(object)) == 0) && is.null(hvg)) {
    warning(paste(
      "No HVGs identified in the object nor provided.",
      "Please run find_hvg_sc() or provide the indices of the HVG",
      "Returning object as is."
    ))
    return(object)
  }

  selected_hvg <- if (!is.null(hvg)) {
    if (.verbose) {
      message(
        paste(
          "HVGs provided.",
          "Will use these ones and set the internal HVG to the provided genes."
        )
      )
    }
    # the method here uses the R 1 indices
    object <- set_hvg(object, hvg)
    hvg
  } else {
    get_hvg(object)
  }

  zeallot::`%<-%`(
    c(pca_factors, pca_loadings, singular_values, scaled),
    rs_sc_pca(
      f_path_gene = get_rust_count_gene_f_path(object),
      no_pcs = no_pcs,
      random_svd = randomised_svd,
      cell_indices = get_cells_to_keep(object),
      gene_indices = selected_hvg,
      seed = seed,
      return_scaled = FALSE,
      verbose = .verbose
    )
  )

  object <- set_pca_factors(object, pca_factors)
  object <- set_pca_loadings(object, pca_loadings)
  object <- set_pca_singular_vals(object, singular_values[1:no_pcs])

  return(object)
}
