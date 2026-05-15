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

### meta cell diffusion coords -------------------------------------------------

#' Calculate diffusion coordinates
#'
#' @description
#' To leverage the quality metrics from Persad, et al., we need the diffusion
#' coordinates to then calculate if a cell is a dense or sparse region of the
#' manifold, its compactness and separation to other meta cells. To do so,
#' generate a diffusion map on the original data based on the approach of
#' SEACells and add the data to the object, see Persad, et al.
#'
#' @param object `MetaCells` class.
#' @param knn_data `SingleCellNearestNeighbour` class. Contains the kNN graph
#' from the original cells.
#' @param n_dcs Integer. Number of diffusion coordinates to use. Defaults to
#' `10L`.
#' @param k_density Integer. The k-th neighbour to use for the density region
#' estimation. Defaults to `150L`.
#' @param seed Integer. Seed for reproducibility
#' @param .verbose Boolean or integer. Controls verbosity and returns run times.
#' `FALSE` -> quiet, `TRUE` or `1L` -> normal verbosity, `2L` -> detailed
#' verbosity.
#'
#' @returns The class with the diffusion map coordinates, density distance and
#' region attached.
#'
#' @export
#'
#' @references
#' Persad, et al. Nat Biotechnol, 2023
calc_diffusion_coordinates <- S7::new_generic(
  name = "calc_diffusion_coordinates",
  dispatch_args = "object",
  fun = function(
    object,
    knn_data,
    n_dcs = 10L,
    k_density = 150L,
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method calc_diffusion_coordinates MetaCells
S7::method(calc_diffusion_coordinates, MetaCells) <- function(
  object,
  knn_data,
  n_dcs = 10L,
  k_density = 150L,
  seed = 42L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, MetaCells))
  checkmate::assertClass(knn_data, "SingleCellNearestNeighbour")
  checkmate::qassert(n_dcs, "I1")
  checkmate::qassert(k_density, "I1")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

  # deal with knn params here
  knn_params <- params_knn_defaults()
  knn_params$k <- 1L

  res <- rs_metacell_density(
    knn_data = knn_data,
    n_dcs = n_dcs,
    k_density = k_density,
    knn_params = knn_params,
    verbose = parse_verbosity(.verbose),
    seed = seed
  )

  density_region <- purrr::map_chr(
    object[[]]$original_cell_idx,
    function(idx) {
      region <- res$regions[idx]
      names(which.max(table(region)))
    }
  )

  object[["density_region"]] <- density_region

  S7::prop(object = object, name = "other_data")[["dcs"]] <- res$dcs
  S7::prop(object = object, name = "other_data")[[
    "density_dist"
  ]] <- res$density_distances
  S7::prop(object = object, name = "other_data")[["regions"]] <- res$regions

  return(object)
}

### manifold metrics -----------------------------------------------------------

#' Calculate manifold metrics
#'
#' @description
#' This function will calculate the compactness and separation of your metacells
#' on the manifold (defined by the diffusion map). You must have run
#' [calc_diffusion_coordinates()] before calling this function. The idea is
#' that compactness indicates how tight the metacell spans the manifold, whereas
#' separation indicates how well the different metacells span the manifold.
#'
#' @returns The class with the compactness and separation scores added.
#'
#' @export
#'
#' @references
#' Persad, et al. Nat Biotechnol, 2023
calc_manifold_metrics <- S7::new_generic(
  name = "calc_manifold_metrics",
  dispatch_args = "object",
  fun = function(
    object
  ) {
    S7::S7_dispatch()
  }
)

#' @method calc_manifold_metrics MetaCells
S7::method(calc_manifold_metrics, MetaCells) <- function(
  object
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, MetaCells))

  dcs <- S7::prop(object, "other_data")[["dcs"]]

  if (is.null(dcs)) {
    warning(
      paste(
        "No diffusion coordinates found.",
        "Please run calc_diffusion_coordinates().",
        "Returning object as is."
      )
    )

    return(object)
  }

  # calculate compactness
  compactness <- rs_metacell_compactness(
    dc = dcs,
    object[[]]$original_cell_idx
  )

  # separation
  separation <- rs_metacell_separation(
    dc = dcs,
    object[[]]$original_cell_idx
  )

  object[["compactness"]] <- compactness
  object[["separation"]] <- separation

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
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

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
  checkmate::qassert(.verbose, c("B1", "I1[0,2]"))

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

  count_list <- mc_counts_to_list(
    object = object,
    gene_indices = selected_hvg,
    assay = "norm"
  )

  zeallot::`%<-%`(
    c(pca_factors, pca_loadings, singular_values),
    rs_mc_pca(
      sparse_data = count_list,
      no_pcs = no_pcs,
      random_svd = randomised_svd,
      seed = seed
    )
  )

  object <- set_pca_factors(object, pca_factors)
  object <- set_pca_loadings(object, pca_loadings)
  object <- set_pca_singular_vals(object, singular_values[1:no_pcs])

  return(object)
}
