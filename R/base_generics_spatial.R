# generics for spatial classes -------------------------------------------------

# generics shared between SpatialSpot, MetaSpot and the SpCache S3 helper.
# spatial-specific S3 generics that mirror the layout of base_generics_sc.R
# but for per-sample state.

## sp_cache setters ------------------------------------------------------------

#' Set the per-sample spatial graph
#'
#' @param x An object with a `SpCache` slot or an `SpCache` itself.
#' @param exp_id String. The experiment identifier.
#' @param graph A `dgCMatrix` (CSR-style sparse matrix) representing the
#' spatial graph for this sample.
#' @param ... Other parameters.
#'
#' @export
#'
#' @keywords internal
set_per_sample_spatial_graph <- function(x, exp_id, graph, ...) {
  UseMethod("set_per_sample_spatial_graph")
}

#' Set the per-sample PCA
#'
#' @param x An object with a `SpCache` slot or an `SpCache` itself.
#' @param exp_id String. The experiment identifier.
#' @param pca Named list with `factors`, `loadings` and `singular_vals`.
#' @param ... Other parameters.
#'
#' @export
#'
#' @keywords internal
set_per_sample_pca <- function(x, exp_id, pca, ...) {
  UseMethod("set_per_sample_pca")
}

#' Set the per-sample Moran's I results
#'
#' @param x An object with a `SpCache` slot or an `SpCache` itself.
#' @param exp_id String. The experiment identifier.
#' @param morans_i `data.table` of Moran's I results for this sample.
#' @param ... Other parameters.
#'
#' @export
#'
#' @keywords internal
set_per_sample_morans_i <- function(x, exp_id, morans_i, ...) {
  UseMethod("set_per_sample_morans_i")
}

#' Set the per-sample SPARK-X results
#'
#' @param x An object with a `SpCache` slot or an `SpCache` itself.
#' @param exp_id String. The experiment identifier.
#' @param sparkx `data.table` of SPARK-X results for this sample.
#' @param ... Other parameters.
#'
#' @export
#'
#' @keywords internal
set_per_sample_sparkx <- function(x, exp_id, sparkx, ...) {
  UseMethod("set_per_sample_sparkx")
}

#' Set the per-sample neighbourhood enrichment results
#'
#' @param x An object with a `SpCache` slot or an `SpCache` itself.
#' @param exp_id String. The experiment identifier.
#' @param label_col String. The obs column name the enrichment was
#' computed against.
#' @param result Named list / `data.table` with the enrichment result.
#' @param ... Other parameters.
#'
#' @export
#'
#' @keywords internal
set_per_sample_nhood_enrichment <- function(
  x,
  exp_id,
  label_col,
  result,
  ...
) {
  UseMethod("set_per_sample_nhood_enrichment")
}

## sp_cache getters ------------------------------------------------------------

#' Get the per-sample spatial graph
#'
#' @param x An object with a `SpCache` slot or an `SpCache` itself.
#' @param exp_id String. The experiment identifier.
#' @param ... Other parameters.
#'
#' @export
get_per_sample_spatial_graph <- function(x, exp_id, ...) {
  UseMethod("get_per_sample_spatial_graph")
}

#' Get the per-sample PCA
#'
#' @param x An object with a `SpCache` slot or an `SpCache` itself.
#' @param exp_id String. The experiment identifier.
#' @param ... Other parameters.
#'
#' @export
get_per_sample_pca <- function(x, exp_id, ...) {
  UseMethod("get_per_sample_pca")
}

#' Get the per-sample Moran's I results
#'
#' @param x An object with a `SpCache` slot or an `SpCache` itself.
#' @param exp_id String. The experiment identifier.
#' @param ... Other parameters.
#'
#' @export
get_per_sample_morans_i <- function(x, exp_id, ...) {
  UseMethod("get_per_sample_morans_i")
}

#' Get the per-sample SPARK-X results
#'
#' @param x An object with a `SpCache` slot or an `SpCache` itself.
#' @param exp_id String. The experiment identifier.
#' @param ... Other parameters.
#'
#' @export
get_per_sample_sparkx <- function(x, exp_id, ...) {
  UseMethod("get_per_sample_sparkx")
}

#' Get the per-sample neighbourhood enrichment results
#'
#' @param x An object with a `SpCache` slot or an `SpCache` itself.
#' @param exp_id String. The experiment identifier.
#' @param label_col String. The obs column name the enrichment was
#' computed against.
#' @param ... Other parameters.
#'
#' @export
get_per_sample_nhood_enrichment <- function(x, exp_id, label_col, ...) {
  UseMethod("get_per_sample_nhood_enrichment")
}

## S7 generics (SpatialSpot / MetaSpot) ----------------------------------------

# Declared here so file load order does not matter. Methods are registered
# in classes_spatial_spot.R and classes_meta_spot.R.

#' Get the samples list from a spatial object
#'
#' @param object A [bixverse::SpatialSpot()] or [bixverse::MetaSpot()].
#'
#' @return Named list of [bixverse::new_spatial_sample()] objects keyed
#' by `exp_id`.
#'
#' @export
get_samples <- S7::new_generic(
  name = "get_samples",
  dispatch_args = "object",
  fun = function(object) {
    S7::S7_dispatch()
  }
)

#' Get a single SpatialSample by exp_id
#'
#' @param object A [bixverse::SpatialSpot()].
#' @param exp_id String. The experiment identifier.
#'
#' @return A [bixverse::new_spatial_sample()].
#'
#' @export
get_sample <- S7::new_generic(
  name = "get_sample",
  dispatch_args = "object",
  fun = function(object, exp_id) {
    S7::S7_dispatch()
  }
)

#' Get the exp_ids present in the samples list
#'
#' @param object A [bixverse::SpatialSpot()] or [bixverse::MetaSpot()].
#'
#' @return Character vector of `exp_id` values.
#'
#' @export
get_sample_ids <- S7::new_generic(
  name = "get_sample_ids",
  dispatch_args = "object",
  fun = function(object) {
    S7::S7_dispatch()
  }
)

#' Get the SpCache from a spatial object
#'
#' @param object A [bixverse::SpatialSpot()].
#'
#' @return The [bixverse::new_sp_cache()] S3 object.
#'
#' @export
get_sp_cache <- S7::new_generic(
  name = "get_sp_cache",
  dispatch_args = "object",
  fun = function(object) {
    S7::S7_dispatch()
  }
)

#' Get spot indices for a given experiment
#'
#' @description
#' Returns 0-based indices (Rust-indexed) for spots belonging to the
#' given experiment, in the ORIGINAL obs space - i.e. positions in the
#' original obs table NOT renumbered after filtering. This matches the
#' indexing convention used by the Rust streaming reader (chunks are
#' keyed by `original_index`).
#'
#' @param object A [bixverse::SpatialSpot()].
#' @param exp_id String. The experiment identifier.
#' @param filtered Boolean. Honour the `to_keep` filter. Defaults to
#' `TRUE`.
#'
#' @return Integer vector of 0-based original obs indices.
#'
#' @export
get_spot_indices_for_exp <- S7::new_generic(
  name = "get_spot_indices_for_exp",
  dispatch_args = "object",
  fun = function(object, exp_id, filtered = TRUE) {
    S7::S7_dispatch()
  }
)

#' Get spatial coordinates for an experiment
#'
#' @param object A [bixverse::SpatialSpot()] or [bixverse::MetaSpot()].
#' @param exp_id String. The experiment identifier.
#' @param filtered Boolean. Honour the `to_keep` filter (SpatialSpot
#' only; ignored for MetaSpot). Defaults to `TRUE`.
#'
#' @return A numeric matrix with two columns `x` and `y`.
#'
#' @export
get_spatial_coords <- S7::new_generic(
  name = "get_spatial_coords",
  dispatch_args = "object",
  fun = function(object, exp_id, filtered = TRUE) {
    S7::S7_dispatch()
  }
)

#' Get the spatial image for an experiment
#'
#' @param object A [bixverse::SpatialSpot()].
#' @param exp_id String. The experiment identifier.
#' @param resolution String. One of `c("lowres", "hires", "fullres")`.
#'
#' @return A numeric array as returned by [png::readPNG()] or
#' [jpeg::readJPEG()].
#'
#' @export
get_image <- S7::new_generic(
  name = "get_image",
  dispatch_args = "object",
  fun = function(
    object,
    exp_id,
    resolution = c("lowres", "hires", "fullres")
  ) {
    S7::S7_dispatch()
  }
)

#' Add a SpatialSample to a spatial object
#'
#' @param object A [bixverse::SpatialSpot()].
#' @param sample A [bixverse::new_spatial_sample()].
#'
#' @return The updated object.
#'
#' @export
add_spatial_sample <- S7::new_generic(
  name = "add_spatial_sample",
  dispatch_args = "object",
  fun = function(object, sample) {
    S7::S7_dispatch()
  }
)

#' Get the per-sample graph from a MetaSpot
#'
#' @param object A [bixverse::MetaSpot()].
#' @param exp_id String. The experiment identifier.
#'
#' @return A `dgCMatrix` adjacency at metaspot resolution.
#'
#' @export
get_per_sample_graph <- S7::new_generic(
  name = "get_per_sample_graph",
  dispatch_args = "object",
  fun = function(object, exp_id) {
    S7::S7_dispatch()
  }
)
