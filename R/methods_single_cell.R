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
#' @param sc_qc_param List. Output of [bixverse::params_sc_min_quality()]. A
#' list with the following elements:
#' \itemize{
#'   \item min_unique_genes - Integer. Minimum number of genes to be detected
#'   in the cell to be included.
#'   \item min_lib_size - Integer. Minimum library size in the cell to be
#'   included.
#'   \item min_cells - Integer. Minimum number of cells a gene needs to be
#'   detected to be included.
#'   \item target_size - Float. Target size to normalise to. Defaults to `1e5`.
#' }
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
    sc_qc_param = params_sc_min_quality(),
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
  sc_qc_param = params_sc_min_quality(),
  .verbose = TRUE
) {
  # checks
  checkmate::assertClass(object, "bixverse::single_cell_exp")
  checkmate::assertFileExists(h5_path)
  assertScMinQC(sc_qc_param)
  checkmate::qassert(.verbose, "B1")

  # rust part
  h5_meta <- get_h5ad_dimensions(h5_path)

  rust_con <- get_sc_rust_ptr(object)

  file_res <- rust_con$h5_to_file(
    cs_type = h5_meta$type,
    h5_path = path.expand(h5_path),
    no_cells = h5_meta$dims["obs"],
    no_genes = h5_meta$dims["var"],
    qc_params = sc_qc_param,
    verbose = .verbose
  )

  rust_con$generate_gene_based_data(
    verbose = .verbose
  )

  # duck db part
  duckdb_con <- get_sc_duckdb(object)
  duckdb_con$populate_obs_from_h5(
    h5_path = h5_path,
    filter = as.integer(file_res$cell_indices + 1)
  )
  duckdb_con$populate_vars_from_h5(
    h5_path = h5_path,
    filter = as.integer(file_res$gene_indices + 1)
  )

  duckdb_con$get_obs_table()

  cell_res_dt <- data.table::setDT(file_res[c("nnz", "lib_size")])

  duckdb_con$add_data_obs(new_data = cell_res_dt)
  cell_map <- duckdb_con$get_obs_index_map()
  gene_map <- duckdb_con$get_var_index_map()

  S7::prop(object, "dims") <- as.integer(rust_con$get_shape())
  object <- set_cell_mapping(x = object, cell_map = cell_map)
  object <- set_gene_mapping(x = object, gene_map = gene_map)

  return(object)
}

### mtx ------------------------------------------------------------------------

#' Load in mtx/plain text files to `single_cell_exp` (nightly!)
#'
#' @description
#' This is a helper function to load in mtx files and corresponding plain text
#' files. It will automatically filter out low quality cells and only keep
#' high quality cells. Under the hood DucKDB and high performance Rust binary
#' files are being used to store the counts.
#'
#' @param object `single_cell_exp` class.
#' @param mtx_path File path to the mtx file you wish to load in.
#' @param obs_path File path
#' @param var_path File path
#' @param sc_qc_param List. Output of [bixverse::params_sc_min_quality()]. A
#' list with the following elements:
#' \itemize{
#'   \item min_unique_genes - Integer. Minimum number of genes to be detected
#'   in the cell to be included.
#'   \item min_lib_size - Integer. Minimum library size in the cell to be
#'   included.
#'   \item min_cells - Integer. Minimum number of cells a gene needs to be
#'   detected to be included.
#'   \item target_size - Float. Target size to normalise to. Defaults to `1e5`.
#' }
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @return It will populate the files on disk and return the class with updated
#' shape information.
#'
#' @export
load_mtx <- S7::new_generic(
  name = "load_mtx",
  dispatch_args = "object",
  fun = function(
    object,
    mtx_path,
    obs_path,
    var_path,
    sc_qc_param = params_sc_min_quality(),
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method load_mtx single_cell_exp
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(load_mtx, single_cell_exp) <- function(
  object,
  mtx_path,
  obs_path,
  var_path,
  sc_qc_param = params_sc_min_quality(),
  .verbose = TRUE
) {
  # checks
  checkmate::assertClass(object, "bixverse::single_cell_exp")
  checkmate::assertFileExists(mtx_path)
  checkmate::assertFileExists(obs_path)
  checkmate::assertFileExists(var_path)
  assertScMinQC(sc_qc_param)
  checkmate::qassert(.verbose, "B1")

  # rust part
  rust_con <- get_sc_rust_ptr(object)

  file_res <- rust_con$mtx_to_file(
    mtx_path = mtx_path,
    qc_params = sc_qc_param,
    verbose = .verbose
  )

  rust_con$generate_gene_based_data(
    verbose = .verbose
  )

  # duckDB part
  duckdb_con <- get_sc_duckdb(object)
  duckdb_con$populate_obs_from_plain_text(
    f_path = obs_path,
    filter = as.integer(file_res$cell_indices + 1)
  )
  duckdb_con$populate_var_from_plain_text(
    f_path = var_path,
    filter = as.integer(file_res$gene_indices + 1)
  )

  duckdb_con$get_obs_table()

  duckdb_con$get_vars_table()

  cell_res_dt <- data.table::setDT(file_res[c("nnz", "lib_size")])

  duckdb_con$add_data_obs(new_data = cell_res_dt)
  cell_map <- duckdb_con$get_obs_index_map()
  gene_map <- duckdb_con$get_var_index_map()

  S7::prop(object, "dims") <- as.integer(rust_con$get_shape())
  object <- set_cell_mapping(x = object, cell_map = cell_map)
  object <- set_gene_mapping(x = object, gene_map = gene_map)

  return(object)
}

## qc --------------------------------------------------------------------------

### gene proportions -----------------------------------------------------------

#' Calculate the proportions of reads for specific gene sets
#'
#' @description
#' This is a helper function that calculates proportions of reads belonging to
#' given gene sets. This can be used for example for the calculation of
#' percentage mitochondrial reads per cell. These will be automatically added
#' to the obs table
#'
#' @param object `single_cell_exp` class.
#' @param gene_set_list A named list with each element containing the gene
#' identifiers of that set. These should be the same as
#' `get_gene_names(object)`!
gene_set_proportions <- S7::new_generic(
  name = "gene_set_proportions",
  dispatch_args = "object",
  fun = function(
    object,
    gene_set_list
  ) {
    S7::S7_dispatch()
  }
)

#' @method gene_set_proportions single_cell_exp
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(gene_set_proportions, single_cell_exp) <- function(
  object,
  gene_set_list
) {
  checkmate::assertClass(object, "bixverse::single_cell_exp")
  checkmate::assertList(gene_set_list, names = "named", types = "character")

  gene_set_list_tidy <- purrr::map(gene_set_list, \(g) {
    get_gene_indices(object, gene_ids = g, rust_index = TRUE)
  })
  names(gene_set_list_tidy) <- names(gene_set_list)

  rs_results <- rs_sc_get_gene_set_perc(
    f_path_cell = file.path(object@dir_data, "counts_cells.bin"),
    gene_set_idx = gene_set_list_tidy
  )

  object[[names(rs_results)]] <- rs_results

  return(object)
}

### hvg ------------------------------------------------------------------------
