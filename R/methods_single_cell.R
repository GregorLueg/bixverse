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
#'
#' @return It will add the columns based on the names in the `gene_set_list` to
#' the obs table.
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
    f_path_cell = get_rust_count_cell_f_path(object),
    gene_set_idx = gene_set_list_tidy
  )

  object[[names(rs_results)]] <- rs_results

  return(object)
}

### hvg ------------------------------------------------------------------------

#' Identify HVGs
#'
#' @description
#' This is a helper function to identify highly variable genes. At the moment
#' the implementation has only
#'
#' @param object `single_cell_exp` class.
#' @param hvg_no Integer. Number of highly variable genes to include. Defaults
#' to `2000L`.
#' @param hvg_params List, see [bixverse::params_sc_hvg()]. This list contains
#' \itemize{
#'   \item method - Which method to use. One of
#'   `c("vst", "meanvarbin", "dispersion")`
#'   \item loess_span - The span for the loess function to standardise the
#'   variance
#'   \item num_bin - Integer. Not yet implemented.
#'   \item bin_method - String. One of `c("equal_width", "equal_freq")`. Not
#'   implemented yet.
#' }
#' @param .verbose Boolean. Controls verbosity and returns run times.
#'
#' @return It will add the columns based on the names in the `gene_set_list` to
#' the obs table.
find_hvg <- S7::new_generic(
  name = "find_hvg",
  dispatch_args = "object",
  fun = function(
    object,
    hvg_no = 2000L,
    hvg_params = params_sc_hvg(),
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method find_hvg single_cell_exp
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(find_hvg, single_cell_exp) <- function(
  object,
  hvg_no = 2000L,
  hvg_params = params_sc_hvg(),
  .verbose = TRUE
) {
  checkmate::assertClass(object, "bixverse::single_cell_exp")
  checkmate::qassert(hvg_no, "I1")
  assertScHvg(hvg_params)
  checkmate::qassert(.verbose, "B1")

  res <- with(
    hvg_params,
    rs_sc_hvg(
      f_path_gene = get_rust_count_gene_f_path(object),
      hvg_method = method,
      cell_indices = get_cells_to_keep(object),
      loess_span = loess_span,
      clip_max = NULL,
      verbose = .verbose
    )
  )

  object <- set_sc_new_var_cols(object = object, data_list = res)

  hvg <- order(res$var_std, decreasing = TRUE)[1:hvg_no]

  object <- set_hvg(object, hvg = hvg)

  return(object)
}

### pca ------------------------------------------------------------------------

#' Run PCA for single cell
#'
#' @description
#' This function will run PCA (option of full SVD and randomised SVD for now)
#' on the detected highly variable genes.
#'
#' @param object `single_cell_exp` class.
#' @param no_pcs Integer. Number of PCs to calculate.
#' @param randomised_svd Boolean. Shall randomised SVD be used. Faster, but
#' less precise.
#' @param seed Integer. Controls reproducibility. Only relevant if
#' `randomised_svd = TRUE`.
#' @param .verbose Boolean. Controls verbosity and returns run times.
#'
#' @return It will add the columns based on the names in the `gene_set_list` to
#' the obs table.
calculate_pca_single_cell <- S7::new_generic(
  name = "calculate_pca_single_cell",
  dispatch_args = "object",
  fun = function(
    object,
    no_pcs,
    randomised_svd = TRUE,
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method calculate_pca_single_cell single_cell_exp
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(calculate_pca_single_cell, single_cell_exp) <- function(
  object,
  no_pcs,
  randomised_svd = TRUE,
  seed = 42L,
  .verbose = TRUE
) {
  checkmate::assertClass(object, "bixverse::single_cell_exp")
  checkmate::qassert(no_pcs, "I1")
  checkmate::qassert(randomised_svd, "B1")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, "B1")

  if (length(get_hvg(object)) == 0) {
    warning(paste(
      "No HVGs identified in this object. Did you run find_hvg()?",
      "Returning object as is."
    ))
    return(object)
  }

  c(pca_factors, pca_loadings) %<-%
    rs_sc_pca(
      f_path_gene = get_rust_count_gene_f_path(object),
      no_pcs = no_pcs,
      random_svd = randomised_svd,
      cell_indices = get_cells_to_keep(object),
      gene_indices = get_hvg(object),
      seed = seed,
      verbose = .verbose
    )

  object <- set_pca_factors(object, pca_factors)
  object <- set_pca_loadings(object, pca_loadings)

  return(object)
}
