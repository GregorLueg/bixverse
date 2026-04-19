# generics that are shared across various S3/S7 classes. the more specific
# generics are found within the given classes_xx.R files.

# shared generics --------------------------------------------------------------

## plotting --------------------------------------------------------------------

#' @title Plot the resolution results.
#'
#' @description
#' Plots the resolution results (if they can be found in the class). The x-axis
#' reflects the different resolutions and the y axis the modularity observed
#' with that resolution.
#'
#' @param object The class, either `RbhGraph` or `BulkCoExp`.
#' @param print_head Boolean. Print the Top5 resolution parameters and their
#' meta data. Only applicable for `BulkCoExp` objects.
#' @param ... Additional arguments to parse to the functions.
#'
#' @return Plots the result, if the results were found in the class. Otherwise,
#' throws a warning and returns NULL.
#'
#' @export
plot_resolution_res <- S7::new_generic(
  name = "plot_resolution_res",
  dispatch_args = "object",
  fun = function(object, print_head = TRUE, ...) {
    S7::S7_dispatch()
  }
)

## meta data -------------------------------------------------------------------

#' @title Replace the meta data
#'
#' @description
#' This function will replace the meta data within the given object
#'
#' @param object The class
#' @param new_metadata data.table. The new meta data you wish to add.
#' @param ... Additional arguments to parse to the functions.
#'
#' @return The object with updated metadata.
#'
#' @export
add_new_metadata <- S7::new_generic(
  name = "add_new_metadata",
  dispatch_args = "object",
  fun = function(object, new_metadata, ...) {
    S7::S7_dispatch()
  }
)
