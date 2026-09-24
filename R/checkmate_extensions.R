# checkmate extensions for various parameter lists provided to the many methods
# in this package

# internal helpers -------------------------------------------------------------

## checks ----------------------------------------------------------------------

### others ---------------------------------------------------------------------

#' Check that files exist
#'
#' @description Checkmate extension for checking if files exist in the
#' directory.
#'
#' @param x String. The directory to check the files for.
#' @param file_names String. Vector of names of the expected files in this
#' directory.
#'
#' @returns `TRUE` if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkFilesExist <- function(x, file_names) {
  res <- purrr::map(file_names, \(file) {
    checkmate::checkFileExists(file.path(x, file))
  })
  res <- purrr::keep(res, \(r) !is.logical(r))
  if (length(res) == 0) {
    return(TRUE)
  }
  res[[1]]
}

#' Assert that files exist
#'
#' @description Checkmate extension for asserting if files exist in the
#' directory.
#'
#' @inheritParams checkFilesExist
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @returns Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertFileExists <- checkmate::makeAssertionFunction(checkFilesExist)

### single cell ----------------------------------------------------------------

#### residuals -----------------------------------------------------------------

#' Check that a clipping range is whole
#'
#' @description Rust falls back to its own default when only one end is given,
#' so half a range silently becomes no range at all.
#'
#' @param x The parameter list.
#' @param label Short label used in the error message.
#'
#' @returns `TRUE` if the check was successful, otherwise an error message.
#'
#' @keywords internal
check_clip_pair <- function(x, label) {
  if (is.null(x[["clip_min"]]) != is.null(x[["clip_max"]])) {
    return(sprintf(
      paste(
        "`clip_min` and `clip_max` in %s must be given together.",
        "Supply both to set a clipping range, or neither to use the default",
        "of +/- sqrt(n_cells)."
      ),
      label
    ))
  }

  if (!is.null(x[["clip_min"]]) && x[["clip_min"]] >= x[["clip_max"]]) {
    return(sprintf("`clip_min` in %s must be smaller than `clip_max`.", label))
  }

  return(TRUE)
}

#' Validate a residual clipping range
#'
#' @description
#' Rust falls back to its own default when only one end is supplied, so half a
#' range silently becomes no range at all. Catch it here instead.
#'
#' @param clip_min Float or `NULL`. Lower bound.
#' @param clip_max Float or `NULL`. Upper bound.
#'
#' @returns Invisibly `TRUE`; called for the error.
#'
#' @keywords internal
assert_clip_range <- function(clip_min, clip_max) {
  checkmate::qassert(clip_min, c("N1", "0"))
  checkmate::qassert(clip_max, c("N1", "0"))

  res <- check_clip_pair(
    list(clip_min = clip_min, clip_max = clip_max),
    "the residual parameters"
  )

  if (!isTRUE(res)) {
    stop(res)
  }

  invisible(TRUE)
}

#### cells in object -----------------------------------------------------------

#' Check that the cell name exists in the object
#'
#' @description Checkmate extension for checking if the provided cell names
#' exist in the object.
#'
#' @param x The `SingleCells` or `SingleCellsSubset` object to check/assert.
#' @param cell_names String. The provided cell names.
#'
#' @returns `TRUE` if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkCellsExist <- function(x, cell_names) {
  res <- checkmate::checkMultiClass(
    x,
    c("bixverse::SingleCells", "bixverse::SingleCellsSubset")
  )
  if (!isTRUE(res)) {
    return(res)
  }
  res <- checkmate::qtest(cell_names, "S+")
  if (!isTRUE(res)) {
    return("The cell names need to be a string vector.")
  }
  if (!all(cell_names %in% get_cell_names(x))) {
    return(paste(
      "Some of the provided cell names do not exist in the object.",
      "Please check."
    ))
  }
  TRUE
}

#' Assert that the cell names exist in the object
#'
#' @description Checkmate extension for asserting if the provided cell names
#' exist in the object.
#'
#' @inheritParams checkCellsExist
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @returns Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertCellsExist <- checkmate::makeAssertionFunction(checkCellsExist)

#### miloR ---------------------------------------------------------------------

#' Check MiloR parameters
#'
#' @description Checkmate extension for checking the MiloR parameters.
#'
#' @param x The list to check/assert
#'
#' @returns `TRUE` if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkScMiloR <- function(x) {
  res <- check_list_shape(
    x,
    c("prop", "k_refine", "refinement_strategy", "index_type")
  )
  if (!isTRUE(res)) {
    return(res)
  }

  # Only a subset of the kNN elements travels with these params, so the rules
  # are applied here rather than through checkKnnParams(), which wants them all.
  res <- apply_qtest_rules(
    x,
    list(
      k_refine = "I1[1,)",
      prop = "N1(0,1)",
      k = "I1[0,)",
      n_trees = "I1[1,)",
      search_budget = c("0", "I1[1,)"),
      delta = "N1[0,1]"
    ),
    label = "MiloR params",
    hint = "k_refine must be an integer >= 1; prop must be in (0, 1)."
  )
  if (!isTRUE(res)) {
    return(res)
  }

  apply_choice_rules(
    x,
    list(
      refinement_strategy = c("approximate", "bruteforce", "index"),
      index_type = c("nndescent", "ivf", "hnsw", "annoy", "exhaustive"),
      knn_method = c("kmknn", "hnsw", "annoy", "nndescent", "ivf", "exhaustive"),
      ann_dist = c("euclidean", "cosine")
    ),
    label = "MiloR params"
  )
}

#' Assert MiloR parameters
#'
#' @description Checkmate extension for asserting the MiloR parameters.
#'
#' @inheritParams checkScMiloR
#'
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @returns Invisibly returns the checked object if the assertion is successful.
#'
#' @keywords internal
assertScMiloR <- checkmate::makeAssertionFunction(checkScMiloR)

#### sc type -------------------------------------------------------------------

#' Check cell marker list
#'
#' @description Checkmate extension for checking a cell marker list as returned
#' by [prepare_cell_markers()].
#'
#' @param x The list to check.
#'
#' @returns `TRUE` if the check was successful, otherwise an error message.
#'
#' @keywords internal
checkCellMarkerList <- function(x) {
  res <- checkmate::checkList(x, names = "named", min.len = 1)
  if (!isTRUE(res)) {
    return(res)
  }

  for (nm in names(x)) {
    entry <- x[[nm]]

    res <- checkmate::checkList(entry, names = "named", len = 3)
    if (!isTRUE(res)) {
      return(sprintf("Entry '%s': %s", nm, res))
    }

    res <- checkmate::checkNames(
      names(entry),
      must.include = c("cell_type", "positive_indices", "negative_indices")
    )
    if (!isTRUE(res)) {
      return(sprintf("Entry '%s': %s", nm, res))
    }

    res <- checkmate::checkString(entry$cell_type, min.chars = 1)
    if (!isTRUE(res)) {
      return(sprintf("Entry '%s' cell_type: %s", nm, res))
    }

    res <- checkmate::checkIntegerish(entry$positive_indices, min.len = 1)
    if (!isTRUE(res)) {
      return(sprintf("Entry '%s' positive_indices: %s", nm, res))
    }

    if (!is.null(entry$negative_indices)) {
      res <- checkmate::checkIntegerish(entry$negative_indices, min.len = 1)
      if (!isTRUE(res)) {
        return(sprintf("Entry '%s' negative_indices: %s", nm, res))
      }
    }
  }

  TRUE
}

#' Assert cell marker list
#'
#' @description Checkmate extension for asserting a cell marker list as returned
#' by [prepare_cell_markers()].
#'
#' @param x The list to assert.
#' @param .var.name Name of the checked object to print in assertions. Defaults
#' to the heuristic implemented in checkmate.
#' @param add Collection to store assertion messages. See
#' [checkmate::makeAssertCollection()].
#'
#' @returns Invisibly returns `x` if the assertion is successful.
#'
#' @keywords internal
assertCellMarkerList <- checkmate::makeAssertionFunction(checkCellMarkerList)
