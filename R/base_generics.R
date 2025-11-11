# shared generics --------------------------------------------------------------

# generics that are shared across various S3/S7 classes. the more specific
# generics are found within the given classes_xx.R files.

## plotting --------------------------------------------------------------------

#' @title Plot the resolution results.
#'
#' @description
#' Plots the resolution results (if they can be found in the class). The x-axis
#' reflects the different resolutions and the y axis the modularity observed
#' with that resolution.
#'
#' @param object The class, either `bixverse::rbh_graph` or `bixverse::bulk_coexp`.
#' @param print_head Boolean. Print the Top5 resolution parameters and their
#' meta data. Only applicable for `bulk_coexp` objects.
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

## single cell -----------------------------------------------------------------

### obs table ------------------------------------------------------------------

#' Getter the obs table
#'
#' @param object `single_cell_exp`, `meta_cells` class.
#' @param indices Optional integer vector. The integer positions of the cells
#' to return.
#' @param cols Optional string vector. The columns from the obs table to return.
#' @param filtered Boolean. Whether to return all cells or filtered to `to_keep`
#' cells. Not relevant for `meta_cells`.
#'
#' @return The obs table
#'
#' @export
get_sc_obs <- S7::new_generic(
  name = "get_sc_obs",
  dispatch_args = "object",
  fun = function(
    object,
    indices = NULL,
    cols = NULL,
    filtered = FALSE
  ) {
    S7::S7_dispatch()
  }
)

### var table ------------------------------------------------------------------

#' Getter the var table
#'
#' @param object `single_cell_exp`, `meta_cells` class.
#' @param indices Optional integer vector. The integer positions of the genes
#' to return.
#' @param cols Optional string vector. The columns from the var table to return.
#'
#' @return The vars table
#'
#' @export
get_sc_var <- S7::new_generic(
  name = "get_sc_var",
  dispatch_args = "object",
  fun = function(
    object,
    indices = NULL,
    cols = NULL
  ) {
    S7::S7_dispatch()
  }
)

### counts ---------------------------------------------------------------------

#' Getter the counts
#'
#' @param object `single_cell_exp`, `meta_cells` class.
#' @param assay String. Which slot to return. One of `c("raw", "norm")`.
#' Defaults to `"raw"`.
#' @param return_format String. One of `c("cell", "gene")`. Return data in
#' cell-centric compressed format (CSR) or gene-centric compressed format (CSC).
#' Defaults to `"cell"`. Not relevant for `meta_cells`.
#' @param cell_indices Optional cell indices.
#' @param gene_indices Optional gene indices.
#' @param use_cells_to_keep Boolean. Shall cells to keep be found in the class,
#' shall the counts be reduced to these. Not relevant for `meta_cells`.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @return The counts table
#'
#' @export
get_sc_counts <- S7::new_generic(
  name = "get_sc_counts",
  dispatch_args = "object",
  fun = function(
    object,
    assay = c("raw", "norm"),
    return_format = c("cell", "gene"),
    cell_indices = NULL,
    gene_indices = NULL,
    use_cells_to_keep = TRUE,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

### knn ------------------------------------------------------------------------

#' Get the KNN matrix
#'
#' @param x An object to get the kNN matrix from.
#'
#' @export
get_knn_mat <- function(x) {
  UseMethod("get_knn_mat")
}

#' Get the KNN distance measures
#'
#' @param x An object to get the kNN distances from.
#'
#' @export
get_knn_dist <- function(x) {
  UseMethod("get_knn_dist")
}
