# single cell methods ----------------------------------------------------------

## i/o -------------------------------------------------------------------------

### h5ad -----------------------------------------------------------------------

#' Load in h5ad to `single_cell_exp` (nightly!)
#'
#' @description
#' This function takes an h5ad file and loads the obs and var data into the
#' DuckDB of the `single_cell_exp` class and the counts into a Rust-binarised
#' format for rapid access. During the reading in of the counts, the log CPM
#' transformation will occur automatically.
#'
#' @param object `single_cell_exp` class.
#' @param h5_path File path to the h5ad object you wish to load in.
#' @param min_genes Integer. Minimum number of genes per cell to be considered
#' in subsequent analysis.
#' @param min_cells Integer. Minimum number of cells per gene to be considered
#' in subsequent analysis.
#' @param target_size Numeric. The target size. 1e6 -> CPM, 1e5 is more commonly
#' used in single cell analysis.
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @return It will populate the files on disk and return the class with updated
#' shape information.
#'
#' @export
load_h5ad <- S7::new_generic(
  name = "load_h5ad",
  dispatch_args = "object",
  fun = function(
    object,
    h5_path,
    min_genes = 200L,
    min_cells = 10L,
    target_size = 1e5,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method load_h5ad single_cell_exp
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(load_h5ad, single_cell_exp) <- function(
  object,
  h5_path,
  min_genes = 200L,
  min_cells = 10L,
  target_size = 1e5,
  .verbose = TRUE
) {
  # checks
  checkmate::assertClass(object, "bixverse::single_cell_exp")
  checkmate::assertFileExists(h5_path)
  checkmate::qassert(min_genes, "I1")
  checkmate::qassert(min_cells, "I1")
  checkmate::qassert(.verbose, "B1")

  # load in the obs and vars
  duckdb_con <- get_sc_duckdb(object)
  duckdb_con$populate_obs_from_h5(h5_path = h5_path, .verbose = .verbose)
  duckdb_con$populate_vars_from_h5(h5_path = h5_path, .verbose = .verbose)

  cell_map <- duckdb_con$get_obs_index_map()
  gene_map <- duckdb_con$get_var_index_map()

  # get meta information from the h5 object
  h5_meta <- get_h5ad_dimensions(h5_path)

  rust_con <- get_sc_rust_ptr(object)

  cell_qc <- rust_con$h5_to_file(
    cs_type = h5_meta$type,
    h5_path = path.expand(h5_path),
    no_cells = h5_meta$dims["obs"],
    no_genes = h5_meta$dims["var"],
    target_size = target_size,
    min_genes = min_genes
  ) %>%
    data.table::setDT()

  gene_qc <- rust_con$generate_gene_based_data(min_cells = min_cells) %>%
    data.table::setDT()

  duckdb_con$add_data_obs(new_data = cell_qc)$add_data_var(new_data = gene_qc)

  S7::prop(object, "dims") = c(h5_meta$dims["obs"], h5_meta$dims["var"])
  S7::prop(object, "index_maps")[["cell_map"]] <- cell_map
  S7::prop(object, "index_maps")[["gene_map"]] <- gene_map

  return(object)
}
