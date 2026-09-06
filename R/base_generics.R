# generics that are shared across various S3/S7 classes. the more specific
# generics are found within the given classes_xx.R files.

# shared generics --------------------------------------------------------------

## plotting --------------------------------------------------------------------

#' Plot the resolution results.
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
#' @returns Plots the result, if the results were found in the class. Otherwise,
#' throws a warning and returns NULL.
#'
#' @export
#'
#' @examples
#' # modularity across the resolutions tested on an RBH graph
#' set.seed(123)
#' modules <- data.table::data.table(
#'   origin = rep(c("set_a", "set_b"), each = 20),
#'   module = rep(c("m1", "m2", "m3", "m4"), each = 10),
#'   gene = unlist(replicate(4, sample(letters, 10), simplify = FALSE))
#' )
#' object <- RbhGraph(
#'   modules,
#'   rbh_type = "set",
#'   dataset_col = "origin",
#'   module_col = "module",
#'   value_col = "gene"
#' )
#' object <- generate_rbh_graph(object, minimum_similarity = 0)
#' object <- find_rbh_communities(object, parallel = FALSE, .verbose = FALSE)
#' plot_resolution_res(object)
plot_resolution_res <- S7::new_generic(
  name = "plot_resolution_res",
  dispatch_args = "object",
  fun = function(object, print_head = TRUE, ...) {
    S7::S7_dispatch()
  }
)

## meta data -------------------------------------------------------------------

#' Replace the meta data
#'
#' @description
#' This function will replace the meta data within the given object
#'
#' @param object The class
#' @param new_metadata data.table. The new meta data you wish to add.
#' @param ... Additional arguments to parse to the functions.
#'
#' @returns The object with updated metadata.
#'
#' @export
#'
#' @examples
#' # swap in a metadata table that carries an extra batch column
#' set.seed(42)
#' counts <- matrix(rpois(60, 20), nrow = 10, ncol = 6)
#' rownames(counts) <- sprintf("gene_%i", 1:10)
#' colnames(counts) <- sprintf("sample_%i", 1:6)
#' meta <- data.table::data.table(
#'   sample_id = colnames(counts),
#'   case_control = rep(c("case", "control"), each = 3)
#' )
#' object <- BulkDge(raw_counts = counts, meta_data = meta)
#' new_meta <- data.table::copy(meta)[, batch := rep(c("b1", "b2"), 3)]
#' object <- add_new_metadata(object, new_metadata = new_meta)
#' head(S7::prop(object, "meta_data"))
add_new_metadata <- S7::new_generic(
  name = "add_new_metadata",
  dispatch_args = "object",
  fun = function(object, new_metadata, ...) {
    S7::S7_dispatch()
  }
)

## scores ----------------------------------------------------------------------

#' Get scores
#'
#' @param x An object to get scores from.
#' @param ... Additional arguments passed to methods.
#'
#' @returns The score matrix or data.table held by the object, depending on the
#' class. Methods exist for `ScTypeResults`, `ScTypeCellResults` and
#' `ScDblFinderRes`.
#'
#' @export
#'
#' @examples
#' # the ScType score matrix of cell types by cluster
#' sc <- demo_single_cells()
#' markers <- data.table::data.table(
#'   cell_type = rep(sprintf("cell_type_%i", 1:3), each = 10),
#'   gene_id = sprintf("gene_%02d", 1:30)
#' )
#' cell_markers <- prepare_cell_markers(obj = sc, marker_df = markers)
#' res <- calc_sc_type_scores(sc, cell_marker_list = cell_markers,
#'                            .verbose = FALSE)
#' dim(get_scores(res))
#'
#' unlink(sc@dir_data, recursive = TRUE, force = TRUE)
get_scores <- function(x, ...) {
  UseMethod("get_scores")
}
