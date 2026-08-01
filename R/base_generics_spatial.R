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
#' @param spot_idx Optional integer vector. The 0-based global spot indices the
#' graph was built over, in graph row order, i.e. what
#' [bixverse::get_spot_indices_for_exp()] returned. A hash of it is stored
#' alongside the graph so the downstream statistics can refuse a graph built
#' over a different set of spots. Without it they fall back to a spot count
#' check, which a `to_keep` change of the same size walks straight past.
#' @param ... Other parameters.
#'
#' @export
#'
#' @keywords internal
set_per_sample_spatial_graph <- function(
  x,
  exp_id,
  graph,
  spot_idx = NULL,
  ...
) {
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

#' Set the per-sample histology image features
#'
#' @param x An object with a `SpCache` slot or an `SpCache` itself.
#' @param exp_id String. The experiment identifier.
#' @param resolution String. The image resolution the features were cut from,
#' one of `c("lowres", "hires", "cytassist", "fullres")`.
#' @param features The image feature result, see
#' [bixverse::new_sp_image_features_res()].
#' @param ... Other parameters.
#'
#' @export
#'
#' @keywords internal
set_per_sample_image_features <- function(
  x,
  exp_id,
  resolution,
  features,
  ...
) {
  UseMethod("set_per_sample_image_features")
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

#' Get the per-sample histology image features
#'
#' @param x An object with a `SpCache` slot or an `SpCache` itself.
#' @param exp_id String. The experiment identifier.
#' @param resolution String. The image resolution the features were cut from,
#' one of `c("lowres", "hires", "cytassist", "fullres")`.
#' @param ... Other parameters.
#'
#' @export
get_per_sample_image_features <- function(x, exp_id, resolution, ...) {
  UseMethod("get_per_sample_image_features")
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
#' @description
#' Returns the spot coordinates in physical units for one experiment.
#'
#' @section Row order:
#' For a [bixverse::SpatialSpot()] the rows come back in exactly the order
#' returned by [bixverse::get_spot_indices_for_exp()] called with the same
#' `exp_id` and `filtered`, which is ascending `cell_idx` (equivalently,
#' ascending original obs position). Row `i` of the coordinate matrix
#' therefore pairs with element `i` of the index vector, so counts pulled
#' with those indices align by position. For a [bixverse::MetaSpot()] the
#' rows are in metaspot order as stored at aggregation time.
#'
#' @param object A [bixverse::SpatialSpot()] or [bixverse::MetaSpot()].
#' @param exp_id String. The experiment identifier.
#' @param filtered Boolean. Honour the `to_keep` filter (SpatialSpot
#' only; ignored for MetaSpot). Defaults to `TRUE`.
#'
#' @return A numeric matrix with two columns `x` and `y`, one row per spot,
#' ordered as described under `Row order`.
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
#' @description
#' Loads a registered slide image off disk. Only PNG and JPEG are readable.
#' `cytassist` is accepted so the path can be requested symmetrically with
#' the other resolutions, but Space Ranger writes that one as a TIFF and
#' TIFF reading is not implemented, so it errors on load. Same for a
#' `fullres` slide handed over as a TIFF.
#'
#' @param object A [bixverse::SpatialSpot()].
#' @param exp_id String. The experiment identifier.
#' @param resolution String. One of `c("lowres", "hires", "cytassist",
#' "fullres")`. Must have been registered in the `image_paths` of the
#' corresponding [bixverse::new_spatial_sample()].
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
    resolution = c("lowres", "hires", "cytassist", "fullres")
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

## S7 generics (spatial analysis) ----------------------------------------------

# Methods are registered in methods_sp_analysis.R against the SpOrSpSubset
# union, so a SpatialSpot and a SpatialSpotSubset behave identically apart from
# how `exp_id` is resolved.

#' Build the spatial neighbourhood graph
#'
#' @description
#' Turns spot coordinates into the neighbour and weight lists every other
#' spatial statistic in this package consumes. The result lands in the
#' `per_sample_spatial_graph` slot of the `SpCache` as a sparse adjacency, one
#' entry per sample.
#'
#' @param object A [bixverse::SpatialSpot()] or
#' [bixverse::SpatialSpotSubset()].
#' @param exp_id Optional string. Restrict to one experiment. `NULL` runs over
#' every registered sample. A `SpatialSpotSubset` only ever has one.
#' @param graph_params List. Output of [bixverse::params_sp_graph()].
#' @param knn_params List or `NULL`. Output of [bixverse::params_sc_knn()],
#' only used by `layout = "knn"`.
#' @param seed Integer. Seed for the ANN index.
#' @param .verbose Boolean or integer. Controls verbosity.
#'
#' @return The object with the graph(s) cached.
#'
#' @export
build_spatial_graph_sp <- S7::new_generic(
  name = "build_spatial_graph_sp",
  dispatch_args = "object",
  fun = function(
    object,
    exp_id = NULL,
    graph_params = params_sp_graph(),
    knn_params = NULL,
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' Detect spatially variable genes with Moran's I
#'
#' @description
#' Per-gene Moran's I against the cached spatial graph, tested under the
#' Cliff-Ord normality null. Needs [bixverse::build_spatial_graph_sp()] to have
#' run first. Results land in the `per_sample_morans_i` slot of the `SpCache`.
#'
#' @param object A [bixverse::SpatialSpot()] or
#' [bixverse::SpatialSpotSubset()].
#' @param exp_id Optional string. Restrict to one experiment. `NULL` runs over
#' every registered sample.
#' @param genes Optional character vector or integer vector. Genes to test.
#' Characters are matched against the var table, integers are taken as
#' 1-indexed positions in it. `NULL` tests every gene.
#' @param svg_params List. Output of [bixverse::params_sp_svg()].
#' @param streaming Boolean. Load the gene chunks in batches of 1000 rather
#' than all at once.
#' @param .verbose Boolean or integer. Controls verbosity.
#'
#' @return The object with the Moran's I results cached.
#'
#' @references
#' Cliff & Ord, Spatial Autocorrelation, 1973
#'
#' @export
morans_i_sp <- S7::new_generic(
  name = "morans_i_sp",
  dispatch_args = "object",
  fun = function(
    object,
    exp_id = NULL,
    genes = NULL,
    svg_params = params_sp_svg(),
    streaming = TRUE,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' Detect spatially variable genes with SPARK-X
#'
#' @description
#' Kernel-based non-parametric test for spatial expression patterns. Works
#' straight off the coordinates and needs no graph. Results land in the
#' `per_sample_sparkx` slot of the `SpCache`.
#'
#' @param object A [bixverse::SpatialSpot()] or
#' [bixverse::SpatialSpotSubset()].
#' @param exp_id Optional string. Restrict to one experiment. `NULL` runs over
#' every registered sample.
#' @param genes Optional character vector or integer vector. Genes to test.
#' Characters are matched against the var table, integers are taken as
#' 1-indexed positions in it. `NULL` tests every gene.
#' @param sparkx_params List. Output of [bixverse::params_sp_sparkx()].
#' @param seed Integer. Seed for landmark selection and bandwidth
#' sub-sampling.
#' @param streaming Boolean. Load the gene chunks in batches of 1000 rather
#' than all at once.
#' @param .verbose Boolean or integer. Controls verbosity.
#'
#' @return The object with the SPARK-X results cached.
#'
#' @references
#' Zhu, Sun & Zhou, Genome Biology, 2021
#'
#' @export
sparkx_sp <- S7::new_generic(
  name = "sparkx_sp",
  dispatch_args = "object",
  fun = function(
    object,
    exp_id = NULL,
    genes = NULL,
    sparkx_params = params_sp_sparkx(),
    seed = 42L,
    streaming = TRUE,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' Neighbourhood enrichment between spot labels
#'
#' @description
#' Asks whether each pair of labels shows up as spatial neighbours more or less
#' often than chance. Only the labels are permuted, so the null keeps the
#' spatial structure of the slide. Needs [bixverse::build_spatial_graph_sp()] to
#' have run first. Results land in the `per_sample_nhood_enrichment` slot of the
#' `SpCache`, keyed by `exp_id` and then by `label_col`.
#'
#' @param object A [bixverse::SpatialSpot()] or
#' [bixverse::SpatialSpotSubset()].
#' @param label_col String. Column in the obs table holding the spot labels.
#' @param exp_id Optional string. Restrict to one experiment. `NULL` runs over
#' every registered sample.
#' @param nhood_params List. Output of [bixverse::params_sp_nhood()].
#' @param seed Integer. Seed for the permutations.
#' @param .verbose Boolean or integer. Controls verbosity.
#'
#' @return The object with the enrichment results cached.
#'
#' @export
nhood_enrichment_sp <- S7::new_generic(
  name = "nhood_enrichment_sp",
  dispatch_args = "object",
  fun = function(
    object,
    label_col,
    exp_id = NULL,
    nhood_params = params_sp_nhood(),
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' Per-spot histology image features
#'
#' @description
#' Cuts one tile per spot out of a registered slide image and reduces it to
#' first-order colour statistics, mean optical density and Haralick texture
#' features. Results land in the `per_sample_image_features` slot of the
#' `SpCache`, keyed by `exp_id` and then by `resolution`.
#'
#' @param object A [bixverse::SpatialSpot()] or
#' [bixverse::SpatialSpotSubset()].
#' @param exp_id Optional string. Restrict to one experiment. `NULL` runs over
#' every registered sample.
#' @param resolution String. Which registered image to cut the tiles from. One
#' of `c("hires", "lowres", "cytassist", "fullres")`.
#' @param image_params List. Output of [bixverse::params_sp_image()].
#' @param .verbose Boolean or integer. Controls verbosity.
#'
#' @return The object with the image features cached.
#'
#' @references
#' Ruifrok & Johnston, Anal Quant Cytol Histol, 2001
#'
#' Haralick, Shanmugam & Dinstein, IEEE Trans Syst Man Cybern, 1973
#'
#' @export
image_features_sp <- S7::new_generic(
  name = "image_features_sp",
  dispatch_args = "object",
  fun = function(
    object,
    exp_id = NULL,
    resolution = c("hires", "lowres", "cytassist", "fullres"),
    image_params = params_sp_image(),
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)
