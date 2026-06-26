# SpatialSample S3 -------------------------------------------------------------

# Per-experiment spatial metadata held inside a named list on a SpatialSpot
# (keyed by exp_id). One SpatialSample per slide/section.

## constructor -----------------------------------------------------------------

#' Construct a SpatialSample
#'
#' @description
#' Helper to build the `SpatialSample` S3 class. Stores per-experiment
#' spatial metadata: technology, number of spots, the obs column names
#' holding coordinates, the parsed scale factors and (relative) paths to
#' the images on disk. Image paths are expected to be RELATIVE to the
#' SpatialSpot `dir_data` so the directory is portable.
#'
#' @param exp_id String. The experiment identifier (matches `exp_id` in
#' the obs table of the parent [bixverse::SpatialSpot()]).
#' @param technology String. One of `c("visium", "xenium")`. Visium HD,
#' MERFISH and Stereo-seq will be added in later milestones.
#' @param n_spots Integer. Number of spots in this experiment.
#' @param coord_cols Named character vector of length 2 with names
#' `c("x", "y")` giving the obs column names holding the spatial
#' coordinates in physical units. Defaults to the Visium full-resolution
#' columns.
#' @param array_cols Named character vector of length 2 with names
#' `c("row", "col")` giving the obs column names with array (grid)
#' coordinates. Visium specific; can be empty for Xenium.
#' @param scale_factors Named list. Parsed `scalefactors_json.json`
#' contents for Visium; named list of relevant scaling info for Xenium.
#' Validated via [bixverse::validate_scale_factors()].
#' @param image_paths Named list of file paths RELATIVE to the parent
#' `dir_data`. Valid names are a subset of `c("lowres", "hires",
#' "fullres")`. Empty list for Xenium.
#'
#' @return A `SpatialSample` S3 object.
#'
#' @export
new_spatial_sample <- function(
  exp_id,
  technology = c("visium", "xenium"),
  n_spots,
  coord_cols = c(x = "pxl_col_in_fullres", y = "pxl_row_in_fullres"),
  array_cols = c(row = "array_row", col = "array_col"),
  scale_factors = list(),
  image_paths = list()
) {
  technology <- match.arg(technology)

  # checks
  checkmate::qassert(exp_id, "S1")
  checkmate::qassert(n_spots, "I1[1,)")
  checkmate::qassert(coord_cols, "S2")
  checkmate::assertNames(names(coord_cols), permutation.of = c("x", "y"))
  checkmate::qassert(array_cols, c("S2", "S0", "0"))
  if (length(array_cols) > 0) {
    checkmate::assertNames(
      names(array_cols),
      permutation.of = c("row", "col")
    )
  }
  checkmate::assertList(scale_factors, names = "named")
  checkmate::assertList(image_paths, names = "named")
  if (length(image_paths) > 0) {
    checkmate::assertSubset(
      names(image_paths),
      choices = c("lowres", "hires", "fullres")
    )
    checkmate::qassert(unlist(image_paths), "S+")
  }

  # technology-specific scale factor validation (no-op for empty list)
  if (length(scale_factors) > 0) {
    validate_scale_factors(scale_factors, technology = technology)
  }

  sample <- list(
    exp_id = exp_id,
    technology = technology,
    n_spots = as.integer(n_spots),
    coord_cols = coord_cols,
    array_cols = array_cols,
    scale_factors = scale_factors,
    image_paths = image_paths
  )

  class(sample) <- "SpatialSample"

  return(sample)
}

## validators ------------------------------------------------------------------

#' Validate a SpatialSample object
#'
#' @description
#' Runs structural checks on a `SpatialSample`. Used by
#' [bixverse::assertSpatialSample()] and during construction of
#' [bixverse::SpatialSpot()].
#'
#' @param x A `SpatialSample` object.
#'
#' @return `TRUE` if the structure is valid, otherwise a character message
#' describing the failure.
#'
#' @export
#'
#' @keywords internal
validate_spatial_sample <- function(x) {
  if (!inherits(x, "SpatialSample")) {
    return("Object is not of class 'SpatialSample'.")
  }
  required <- c(
    "exp_id",
    "technology",
    "n_spots",
    "coord_cols",
    "array_cols",
    "scale_factors",
    "image_paths"
  )
  missing_fields <- setdiff(required, names(x))
  if (length(missing_fields) > 0) {
    return(sprintf(
      "SpatialSample missing field(s): %s",
      paste(missing_fields, collapse = ", ")
    ))
  }
  if (!is.character(x$exp_id) || length(x$exp_id) != 1L) {
    return("'exp_id' must be a single string.")
  }
  if (!x$technology %in% c("visium", "xenium")) {
    return("'technology' must be one of c('visium', 'xenium').")
  }
  if (!is.integer(x$n_spots) || x$n_spots < 1L) {
    return("'n_spots' must be a positive integer.")
  }
  if (
    !is.character(x$coord_cols) ||
      length(x$coord_cols) != 2L ||
      !setequal(names(x$coord_cols), c("x", "y"))
  ) {
    return("'coord_cols' must be a length-2 named character with names x, y.")
  }
  return(TRUE)
}

#' Validate spatial scale factors for a given technology
#'
#' @description
#' Enforces the expected keys in a `scale_factors` list. For Visium this
#' matches the contents of `scalefactors_json.json`. For Xenium a small
#' set of pixel-size related keys is required.
#'
#' @param scale_factors Named list. Parsed scale factors.
#' @param technology String. One of `c("visium", "xenium")`.
#'
#' @return Invisibly `TRUE` if the validation passes, otherwise throws an
#' informative error.
#'
#' @export
validate_scale_factors <- function(
  scale_factors,
  technology = c("visium", "xenium")
) {
  technology <- match.arg(technology)
  checkmate::assertList(scale_factors, names = "named")

  required <- switch(
    technology,
    visium = c(
      "spot_diameter_fullres",
      "tissue_hires_scalef",
      "tissue_lowres_scalef",
      "fiducial_diameter_fullres"
    ),
    xenium = c("pixel_size")
  )

  missing_keys <- setdiff(required, names(scale_factors))
  if (length(missing_keys) > 0) {
    stop(sprintf(
      "scale_factors for technology '%s' missing key(s): %s",
      technology,
      paste(missing_keys, collapse = ", ")
    ))
  }

  for (k in required) {
    val <- scale_factors[[k]]
    if (!is.numeric(val) || length(val) != 1L || !is.finite(val)) {
      stop(sprintf(
        "scale_factors[['%s']] must be a single finite numeric value.",
        k
      ))
    }
  }

  invisible(TRUE)
}

## checkmate extension ---------------------------------------------------------

#' Check that an object is a valid SpatialSample
#'
#' @description Checkmate extension for [bixverse::SpatialSample()].
#'
#' @param x The object to check.
#'
#' @return `TRUE` if the check was successful, otherwise an error message
#' string.
#'
#' @keywords internal
checkSpatialSample <- function(x) {
  res <- validate_spatial_sample(x)
  if (isTRUE(res)) TRUE else res
}

#' Assert that an object is a valid SpatialSample
#'
#' @description Checkmate assertion for [bixverse::SpatialSample()].
#'
#' @inheritParams checkSpatialSample
#'
#' @param .var.name Name of the checked object to print in assertions.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @return Invisibly returns the checked object if the assertion is
#' successful.
#'
#' @keywords internal
assertSpatialSample <- checkmate::makeAssertionFunction(checkSpatialSample)

## primitives ------------------------------------------------------------------

#' Print method for SpatialSample objects
#'
#' @param x A `SpatialSample` object.
#' @param ... Unused.
#'
#' @return Invisibly returns `x`.
#'
#' @export
print.SpatialSample <- function(x, ...) {
  cat(
    sprintf("SpatialSample (exp_id = '%s')\n", x$exp_id),
    sprintf("  Technology: %s\n", x$technology),
    sprintf("  No spots: %i\n", x$n_spots),
    sprintf(
      "  Coord cols: x = '%s', y = '%s'\n",
      x$coord_cols[["x"]],
      x$coord_cols[["y"]]
    ),
    sprintf(
      "  Image(s): %s\n",
      if (length(x$image_paths) == 0) {
        "none"
      } else {
        paste(names(x$image_paths), collapse = ", ")
      }
    ),
    sep = ""
  )
  invisible(x)
}
