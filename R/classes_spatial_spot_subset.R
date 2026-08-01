# spatial spot subset class ----------------------------------------------------

## s7 object -------------------------------------------------------------------

#' @title bixverse spatial spot subset class
#'
#' @description
#' Subset view onto a [bixverse::SpatialSpot()] object, restricted to the spots
#' of a single experiment. Mirrors [bixverse::SingleCellsSubset()] and inherits
#' from it, so every single-cell-style getter, setter and method works
#' unchanged. The Rust count connection is shared with the parent (no data
#' copy); `obs_table` and `var_table` are held in memory; `sc_map` is rebuilt to
#' point only at the subset spots but stays in the original index space so Rust
#' calls remain valid without further translation.
#'
#' The grouping column is always `exp_id`, so the constructor takes the
#' experiment identifier rather than a generic column/value pair.
#'
#' @section Properties:
#' \describe{
#'   \item{count_connection}{Shared Rust pointer to the on-disk counts.}
#'   \item{dir_data}{Directory holding the binary count files.}
#'   \item{obs_table}{Subset obs, sorted by `cell_idx`. `cell_idx` keeps the
#'   original 1-indexed position in the parent.}
#'   \item{var_table}{Variable/feature table (unchanged from parent).}
#'   \item{grouping_column}{Always `"exp_id"`.}
#'   \item{group}{The `exp_id` this subset represents.}
#'   \item{sc_cache}{Fresh `ScCache` for subset-specific PCA, kNN, sNN,
#'   embeddings.}
#'   \item{sc_map}{`ScMap` restricted to the subset spots. `cell_mapping`
#'   stays 1-indexed and `cells_to_keep_idx` stays 0-indexed, both in the
#'   original parent index space.}
#'   \item{subset_to_original}{Integer vector. 1-indexed original spot
#'   positions, ascending, in subset row order.}
#'   \item{dims}{`c(n_spots_subset, n_genes)`.}
#'   \item{sample}{The [bixverse::new_spatial_sample()] of this experiment.}
#'   \item{sp_cache}{Fresh [bixverse::new_sp_cache()] for the spatial results.}
#'   \item{coords}{Numeric matrix with columns `x` and `y` in full-resolution
#'   pixels, one row per subset spot, aligned with `subset_to_original`.}
#' }
#'
#' @section Row order:
#' `obs_table`, `coords` and `subset_to_original` are all in ascending
#' `cell_idx` order, which is what [bixverse::get_spot_indices_for_exp()]
#' returns. Row `i` of `coords` therefore belongs to spot
#' `subset_to_original[i]` of the parent, and everything handed to Rust lines up
#' by position without a translation step.
#'
#' @param sp_object A [bixverse::SpatialSpot()].
#' @param exp_id String. The experiment identifier to subset to.
#'
#' @return A `SpatialSpotSubset` object.
#'
#' @export
SpatialSpotSubset <- S7::new_class(
  name = "SpatialSpotSubset",
  parent = SingleCellsSubset,
  properties = list(
    count_connection = S7::class_any,
    dir_data = S7::class_character,
    obs_table = S7::class_data.frame,
    var_table = S7::class_data.frame,
    grouping_column = S7::class_character,
    group = S7::class_character,
    sc_cache = S7::class_any,
    sc_map = S7::class_any,
    subset_to_original = S7::class_integer,
    dims = S7::class_integer,
    sample = S7::class_any,
    sp_cache = S7::class_any,
    coords = S7::class_any
  ),
  constructor = function(sp_object, exp_id) {
    checkmate::assertClass(sp_object, "bixverse::SpatialSpot")
    checkmate::qassert(exp_id, "S1")

    dir_data <- S7::prop(sp_object, "dir_data")
    sc_map <- S7::prop(sp_object, "sc_map")
    count_connection <- sp_object@count_connection
    sample <- get_sample(sp_object, exp_id)

    obs_table <- get_sc_obs(sp_object, filtered = TRUE)
    checkmate::assertNames(x = colnames(obs_table), must.include = "exp_id")
    obs_table <- obs_table[get("exp_id") == exp_id, ]

    if (nrow(obs_table) == 0L) {
      stop(sprintf("No spots found for exp_id '%s'.", exp_id))
    }

    # the DuckDB query behind `get_sc_obs()` has no ORDER BY, but the
    # coordinates come back in ascending `cell_idx`. Sort here or the two stop
    # lining up and nothing downstream notices.
    data.table::setorderv(obs_table, "cell_idx")

    var_table <- get_sc_var(sp_object)

    subset_to_original <- as.integer(obs_table$cell_idx)

    spot_indices <- get_spot_indices_for_exp(
      sp_object,
      exp_id = exp_id,
      filtered = TRUE
    )
    if (!identical(as.integer(spot_indices), subset_to_original - 1L)) {
      stop(sprintf(
        "The obs rows of exp_id '%s' disagree with its spot indices.",
        exp_id
      ))
    }

    coords <- get_spatial_coords(sp_object, exp_id = exp_id, filtered = TRUE)
    if (nrow(coords) != nrow(obs_table)) {
      stop(sprintf(
        "Got %i coordinate rows for %i spots of exp_id '%s'.",
        nrow(coords),
        nrow(obs_table),
        exp_id
      ))
    }

    if (!is.null(sc_map$cells_to_keep_idx)) {
      # cells_to_keep_idx is 0-indexed; compare to 0-indexed subset positions
      sc_map$cells_to_keep_idx <- sc_map$cells_to_keep_idx[
        sc_map$cells_to_keep_idx %in% (subset_to_original - 1L)
      ]
    }
    # HVG is subset-specific
    sc_map$hvg_gene_indices <- NULL

    S7::new_object(
      S7::S7_object(),
      count_connection = count_connection,
      dir_data = dir_data,
      obs_table = obs_table,
      var_table = var_table,
      grouping_column = "exp_id",
      group = exp_id,
      sc_cache = new_sc_cache(),
      sc_map = sc_map,
      subset_to_original = subset_to_original,
      dims = c(nrow(obs_table), nrow(var_table)),
      sample = sample,
      sp_cache = new_sp_cache(),
      coords = coords
    )
  }
)

## primitives ------------------------------------------------------------------

#' @name print.SpatialSpotSubset
#'
#' @title print Method for SpatialSpotSubset object
#'
#' @param x An object of class `SpatialSpotSubset`.
#' @param ... Additional arguments (currently not used).
#'
#' @returns Invisibly returns `x`.
#'
#' @method print SpatialSpotSubset
#'
#' @keywords internal
S7::method(print, SpatialSpotSubset) <- function(x, ...) {
  checkmate::assertTRUE(S7::S7_inherits(x, SpatialSpotSubset))

  dims <- S7::prop(x, "dims")
  exp_id <- S7::prop(x, "group")
  sample <- S7::prop(x, "sample")
  sp_cache <- S7::prop(x, "sp_cache")
  sc_cache <- S7::prop(x, "sc_cache")

  nhood_cols <- names(sp_cache$per_sample_nhood_enrichment[[exp_id]])
  image_res <- names(sp_cache$per_sample_image_features[[exp_id]])

  cat(
    "Spatial transcriptomics experiment (subset).\n",
    sprintf("  Experiment: %s (%s)\n", exp_id, sample$technology),
    sprintf("  No spots: %i\n", dims[1]),
    sprintf("  No genes: %i\n", dims[2]),
    sprintf(
      "  PCA calculated: %s\n",
      !is.null(sc_cache[["pca_factors"]])
    ),
    sprintf(
      "  Spatial graph: %s\n",
      !is.null(sp_cache$per_sample_spatial_graph[[exp_id]])
    ),
    sprintf(
      "  Moran's I: %s\n",
      !is.null(sp_cache$per_sample_morans_i[[exp_id]])
    ),
    sprintf(
      "  SPARK-X: %s\n",
      !is.null(sp_cache$per_sample_sparkx[[exp_id]])
    ),
    sprintf(
      "  Nhood enrichment: %s\n",
      if (length(nhood_cols) > 0) toString(nhood_cols) else "none"
    ),
    sprintf(
      "  Image features: %s\n",
      if (length(image_res) > 0) toString(image_res) else "none"
    ),
    sep = ""
  )
  invisible(x)
}

## getters ---------------------------------------------------------------------

#' Assert an exp_id is the one a subset holds
#'
#' @description
#' A subset carries exactly one experiment, so anything else is a caller
#' mistake rather than a lookup miss.
#'
#' @param object A [bixverse::SpatialSpotSubset()].
#' @param exp_id String. The experiment identifier to check.
#'
#' @return Invisibly `TRUE`, otherwise it errors.
#'
#' @keywords internal
.sp_subset_assert_exp_id <- function(object, exp_id) {
  checkmate::qassert(exp_id, "S1")
  if (!identical(exp_id, S7::prop(object, "group"))) {
    stop(sprintf(
      "The subset holds exp_id '%s', not '%s'.",
      S7::prop(object, "group"),
      exp_id
    ))
  }
  invisible(TRUE)
}

### samples --------------------------------------------------------------------

#' @method get_samples SpatialSpotSubset
#'
#' @export
S7::method(get_samples, SpatialSpotSubset) <- function(object) {
  checkmate::assertTRUE(S7::S7_inherits(object, SpatialSpotSubset))
  samples <- list(S7::prop(object, "sample"))
  names(samples) <- S7::prop(object, "group")
  return(samples)
}

#' @method get_sample SpatialSpotSubset
#'
#' @export
S7::method(get_sample, SpatialSpotSubset) <- function(object, exp_id) {
  checkmate::assertTRUE(S7::S7_inherits(object, SpatialSpotSubset))
  .sp_subset_assert_exp_id(object, exp_id)

  return(S7::prop(object, "sample"))
}

#' @method get_sample_ids SpatialSpotSubset
#'
#' @export
S7::method(get_sample_ids, SpatialSpotSubset) <- function(object) {
  checkmate::assertTRUE(S7::S7_inherits(object, SpatialSpotSubset))
  return(S7::prop(object, "group"))
}

### sp_cache -------------------------------------------------------------------

#' @method get_sp_cache SpatialSpotSubset
#'
#' @export
S7::method(get_sp_cache, SpatialSpotSubset) <- function(object) {
  checkmate::assertTRUE(S7::S7_inherits(object, SpatialSpotSubset))
  return(S7::prop(object, "sp_cache"))
}

### spot indices ---------------------------------------------------------------

#' @method get_spot_indices_for_exp SpatialSpotSubset
#'
#' @export
S7::method(get_spot_indices_for_exp, SpatialSpotSubset) <- function(
  object,
  exp_id,
  filtered = TRUE
) {
  checkmate::assertTRUE(S7::S7_inherits(object, SpatialSpotSubset))
  checkmate::qassert(filtered, "B1")
  .sp_subset_assert_exp_id(object, exp_id)

  # `filtered` kept for parent API parity; the subset is always filtered.
  return(S7::prop(object, "subset_to_original") - 1L)
}

### spatial coords -------------------------------------------------------------

#' @method get_spatial_coords SpatialSpotSubset
#'
#' @export
S7::method(get_spatial_coords, SpatialSpotSubset) <- function(
  object,
  exp_id,
  filtered = TRUE
) {
  checkmate::assertTRUE(S7::S7_inherits(object, SpatialSpotSubset))
  checkmate::qassert(filtered, "B1")
  .sp_subset_assert_exp_id(object, exp_id)

  return(S7::prop(object, "coords"))
}

### images ---------------------------------------------------------------------

#' @method get_image SpatialSpotSubset
#'
#' @export
S7::method(get_image, SpatialSpotSubset) <- function(
  object,
  exp_id,
  resolution = c("lowres", "hires", "cytassist", "fullres")
) {
  resolution <- match.arg(resolution)

  checkmate::assertTRUE(S7::S7_inherits(object, SpatialSpotSubset))
  checkmate::assertChoice(
    resolution,
    c("lowres", "hires", "cytassist", "fullres")
  )

  sample <- get_sample(object, exp_id)
  rel_path <- sample$image_paths[[resolution]]
  if (is.null(rel_path)) {
    stop(sprintf(
      "No image at resolution '%s' registered for exp_id '%s'.",
      resolution,
      exp_id
    ))
  }
  return(.read_spatial_image(
    file.path(S7::prop(object, "dir_data"), rel_path)
  ))
}

## sp_cache forwarding (SpatialSpotSubset -> SpCache) --------------------------

# Same story as on `SpatialSpot`: these are S3 generics targeting an S7 class,
# so they have to go through `S7::method()`.

#' @rdname set_per_sample_spatial_graph
#'
#' @export
S7::method(set_per_sample_spatial_graph, SpatialSpotSubset) <- function(
  x,
  exp_id,
  graph,
  ...
) {
  checkmate::assertTRUE(S7::S7_inherits(x, SpatialSpotSubset))
  S7::prop(x, "sp_cache") <- set_per_sample_spatial_graph(
    x = S7::prop(x, "sp_cache"),
    exp_id = exp_id,
    graph = graph
  )
  return(x)
}

#' @rdname set_per_sample_pca
#'
#' @export
S7::method(set_per_sample_pca, SpatialSpotSubset) <- function(
  x,
  exp_id,
  pca,
  ...
) {
  checkmate::assertTRUE(S7::S7_inherits(x, SpatialSpotSubset))
  S7::prop(x, "sp_cache") <- set_per_sample_pca(
    x = S7::prop(x, "sp_cache"),
    exp_id = exp_id,
    pca = pca
  )
  return(x)
}

#' @rdname set_per_sample_morans_i
#'
#' @export
S7::method(set_per_sample_morans_i, SpatialSpotSubset) <- function(
  x,
  exp_id,
  morans_i,
  ...
) {
  checkmate::assertTRUE(S7::S7_inherits(x, SpatialSpotSubset))
  S7::prop(x, "sp_cache") <- set_per_sample_morans_i(
    x = S7::prop(x, "sp_cache"),
    exp_id = exp_id,
    morans_i = morans_i
  )
  return(x)
}

#' @rdname set_per_sample_sparkx
#'
#' @export
S7::method(set_per_sample_sparkx, SpatialSpotSubset) <- function(
  x,
  exp_id,
  sparkx,
  ...
) {
  checkmate::assertTRUE(S7::S7_inherits(x, SpatialSpotSubset))
  S7::prop(x, "sp_cache") <- set_per_sample_sparkx(
    x = S7::prop(x, "sp_cache"),
    exp_id = exp_id,
    sparkx = sparkx
  )
  return(x)
}

#' @rdname set_per_sample_nhood_enrichment
#'
#' @export
S7::method(set_per_sample_nhood_enrichment, SpatialSpotSubset) <- function(
  x,
  exp_id,
  label_col,
  result,
  ...
) {
  checkmate::assertTRUE(S7::S7_inherits(x, SpatialSpotSubset))
  S7::prop(x, "sp_cache") <- set_per_sample_nhood_enrichment(
    x = S7::prop(x, "sp_cache"),
    exp_id = exp_id,
    label_col = label_col,
    result = result
  )
  return(x)
}

#' @rdname set_per_sample_image_features
#'
#' @export
S7::method(set_per_sample_image_features, SpatialSpotSubset) <- function(
  x,
  exp_id,
  resolution,
  features,
  ...
) {
  checkmate::assertTRUE(S7::S7_inherits(x, SpatialSpotSubset))
  S7::prop(x, "sp_cache") <- set_per_sample_image_features(
    x = S7::prop(x, "sp_cache"),
    exp_id = exp_id,
    resolution = resolution,
    features = features
  )
  return(x)
}

#' @rdname get_per_sample_spatial_graph
#'
#' @export
S7::method(get_per_sample_spatial_graph, SpatialSpotSubset) <- function(
  x,
  exp_id,
  ...
) {
  checkmate::assertTRUE(S7::S7_inherits(x, SpatialSpotSubset))
  return(get_per_sample_spatial_graph(
    x = S7::prop(x, "sp_cache"),
    exp_id = exp_id
  ))
}

#' @rdname get_per_sample_pca
#'
#' @export
S7::method(get_per_sample_pca, SpatialSpotSubset) <- function(x, exp_id, ...) {
  checkmate::assertTRUE(S7::S7_inherits(x, SpatialSpotSubset))
  return(get_per_sample_pca(
    x = S7::prop(x, "sp_cache"),
    exp_id = exp_id
  ))
}

#' @rdname get_per_sample_morans_i
#'
#' @export
S7::method(get_per_sample_morans_i, SpatialSpotSubset) <- function(
  x,
  exp_id,
  ...
) {
  checkmate::assertTRUE(S7::S7_inherits(x, SpatialSpotSubset))
  return(get_per_sample_morans_i(
    x = S7::prop(x, "sp_cache"),
    exp_id = exp_id
  ))
}

#' @rdname get_per_sample_sparkx
#'
#' @export
S7::method(get_per_sample_sparkx, SpatialSpotSubset) <- function(
  x,
  exp_id,
  ...
) {
  checkmate::assertTRUE(S7::S7_inherits(x, SpatialSpotSubset))
  return(get_per_sample_sparkx(
    x = S7::prop(x, "sp_cache"),
    exp_id = exp_id
  ))
}

#' @rdname get_per_sample_nhood_enrichment
#'
#' @export
S7::method(get_per_sample_nhood_enrichment, SpatialSpotSubset) <- function(
  x,
  exp_id,
  label_col,
  ...
) {
  checkmate::assertTRUE(S7::S7_inherits(x, SpatialSpotSubset))
  return(get_per_sample_nhood_enrichment(
    x = S7::prop(x, "sp_cache"),
    exp_id = exp_id,
    label_col = label_col
  ))
}

#' @rdname get_per_sample_image_features
#'
#' @export
S7::method(get_per_sample_image_features, SpatialSpotSubset) <- function(
  x,
  exp_id,
  resolution,
  ...
) {
  checkmate::assertTRUE(S7::S7_inherits(x, SpatialSpotSubset))
  return(get_per_sample_image_features(
    x = S7::prop(x, "sp_cache"),
    exp_id = exp_id,
    resolution = resolution
  ))
}
