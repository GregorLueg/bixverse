# single cell i/o --------------------------------------------------------------

## seurat ----------------------------------------------------------------------

#' Load in Seurat to `SingleCells`
#'
#' @description
#' This function takes a Seurat file and generates `SingleCells` class
#' from it.
#'
#' @param object `SingleCells` class.
#' @param seurat `Seurat` class you want to transform.
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
#' @param streaming Boolean. Shall the data be streamed during the conversion
#' of CSR to CSC. Defaults to `TRUE` and should be used for larger data sets.
#' @param batch_size Integer. If `streaming = TRUE`, how many cells to process
#' in one batch. Defaults to `1000L`.
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @return It will populate the files on disk and return the class with updated
#' shape information.
#'
#' @export
load_seurat <- S7::new_generic(
  name = "load_h5ad",
  dispatch_args = "object",
  fun = function(
    object,
    seurat,
    sc_qc_param = params_sc_min_quality(),
    batch_size = 1000L,
    streaming = TRUE,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method load_seurat SingleCells
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(load_seurat, SingleCells) <- function(
  object,
  seurat,
  sc_qc_param = params_sc_min_quality(),
  batch_size = 1000L,
  streaming = TRUE,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::assertClass(seurat, "Seurat")
  assertScMinQC(sc_qc_param)
  checkmate::qassert(batch_size, "I1")
  checkmate::qassert(streaming, "B1")
  checkmate::qassert(.verbose, "B1")

  counts <- get_seurat_counts_to_list(seurat)
  obs_dt <- data.table::as.data.table(
    seurat@meta.data,
    keep.rownames = "barcode"
  )

  if ("cell_id" %in% names(obs_dt)) {
    obs_dt[, cell_id := NULL]
  }

  var_dt <- data.table::data.table(gene_id = rownames(seurat))

  rust_con <- get_sc_rust_ptr(object)

  file_res <- rust_con$r_data_to_file(
    r_data = counts,
    qc_params = sc_qc_param,
    verbose = .verbose
  )

  if (streaming) {
    rust_con$generate_gene_based_data_streaming(
      batch_size = batch_size,
      verbose = .verbose
    )
  } else {
    rust_con$generate_gene_based_data(
      verbose = .verbose
    )
  }

  gene_nnz <- rust_con$get_nnz_genes(gene_indices = NULL)
  gene_nnz_dt <- data.table::data.table(no_cells_exp = gene_nnz)

  duckdb_con <- get_sc_duckdb(object)

  duckdb_con$populate_obs_from_data.table(
    obs_dt = obs_dt,
    filter = as.integer(file_res$cell_indices + 1)
  )

  duckdb_con$populate_var_from_data.table(
    var_dt = var_dt,
    filter = as.integer(file_res$gene_indices + 1)
  )

  cell_res_dt <- data.table::setDT(file_res[c("nnz", "lib_size")])

  duckdb_con$add_data_obs(new_data = cell_res_dt)
  duckdb_con$add_data_var(new_data = gene_nnz_dt)
  duckdb_con$set_to_keep_column()
  cell_map <- duckdb_con$get_obs_index_map()
  gene_map <- duckdb_con$get_var_index_map()

  S7::prop(object, "dims") <- as.integer(rust_con$get_shape())
  object <- set_cell_mapping(x = object, cell_map = cell_map)
  object <- set_gene_mapping(x = object, gene_map = gene_map)

  return(object)
}

## direct r --------------------------------------------------------------------

#' Load in data directly from R objects.
#'
#' @description
#' This function loads in data directly from R objects.
#'
#' @param object `SingleCells` class.
#' @param counts Sparse matrix. The cells represent the rows, the genes the
#' indices. Needs to be `"dgRMatrix"`.
#' @param obs data.table. The data.table representing the observations, i.e.,
#' cell information.
#' @param var data.table. The data.table representing the features, i.e., the
#' feature information.
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
#' @param streaming Boolean. Shall the data be streamed during the conversion
#' of CSR to CSC. Defaults to `TRUE` and should be used for larger data sets.
#' @param batch_size Integer. If `streaming = TRUE`, how many cells to process
#' in one batch. Defaults to `1000L`.
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @return It will populate the files on disk and return the class with updated
#' shape information.
#'
#' @export
load_r_data <- S7::new_generic(
  name = "load_h5ad",
  dispatch_args = "object",
  fun = function(
    object,
    counts,
    obs,
    var,
    sc_qc_param = params_sc_min_quality(),
    batch_size = 1000L,
    streaming = TRUE,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method load_r_data SingleCells
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(load_r_data, SingleCells) <- function(
  object,
  counts,
  obs,
  var,
  sc_qc_param = params_sc_min_quality(),
  batch_size = 1000L,
  streaming = TRUE,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::assertClass(counts, "dgRMatrix")
  no_cells <- nrow(counts)
  no_genes <- ncol(counts)
  checkmate::assertDataTable(obs, nrows = no_cells)
  checkmate::assertDataTable(var, nrows = no_genes)
  checkmate::qassert(.verbose, "B1")

  if (.verbose) {
    message("Writing counts to disk.")
  }
  # body
  counts <- sparse_mat_to_list(counts)

  rust_con <- get_sc_rust_ptr(object)

  file_res <- rust_con$r_data_to_file(
    r_data = counts,
    qc_params = sc_qc_param,
    verbose = .verbose
  )

  if (.verbose) {
    message("Generating gene-based data.")
  }
  if (streaming) {
    rust_con$generate_gene_based_data_streaming(
      batch_size = batch_size,
      verbose = .verbose
    )
  } else {
    rust_con$generate_gene_based_data(
      verbose = .verbose
    )
  }

  gene_nnz <- rust_con$get_nnz_genes(gene_indices = NULL)
  gene_nnz_dt <- data.table::data.table(no_cells_exp = gene_nnz)

  if (.verbose) {
    message("Writing to the DuckDB.")
  }
  duckdb_con <- get_sc_duckdb(object)

  duckdb_con$populate_obs_from_data.table(
    obs_dt = obs,
    filter = as.integer(file_res$cell_indices + 1)
  )

  duckdb_con$populate_var_from_data.table(
    var_dt = var,
    filter = as.integer(file_res$gene_indices + 1)
  )

  cell_res_dt <- data.table::setDT(file_res[c("nnz", "lib_size")])

  if (.verbose) {
    message("Setting internal mapping.")
  }
  duckdb_con$add_data_obs(new_data = cell_res_dt)
  duckdb_con$add_data_var(new_data = gene_nnz_dt)
  duckdb_con$set_to_keep_column()
  cell_map <- duckdb_con$get_obs_index_map()
  gene_map <- duckdb_con$get_var_index_map()

  S7::prop(object, "dims") <- as.integer(rust_con$get_shape())
  object <- set_cell_mapping(x = object, cell_map = cell_map)
  object <- set_gene_mapping(x = object, gene_map = gene_map)

  return(object)
}

### h5ad -----------------------------------------------------------------------

#### fast ----------------------------------------------------------------------

#' Load in h5ad to `SingleCells`
#'
#' @description
#' This function takes an h5ad file and loads the obs and var data into the
#' DuckDB of the `SingleCells` class and the counts into a Rust-binarised
#' format for rapid access. During the reading in of the counts, the log CPM
#' transformation will occur automatically.
#'
#' @param object `SingleCells` class.
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
#' @param cell_id_col Optional string. If a specific column in the h5ad
#' obs data is representing the cell identifiers, you can specify it here.
#' @param streaming Boolean. Shall the data be streamed during the conversion
#' of CSR to CSC. Defaults to `TRUE` and should be used for larger data sets.
#' @param batch_size Integer. If `streaming = TRUE`, how many cells to process
#' in one batch. Defaults to `1000L`.
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
    streaming = TRUE,
    cell_id_col = NULL,
    batch_size = 1000L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method load_h5ad SingleCells
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(load_h5ad, SingleCells) <- function(
  object,
  h5_path,
  sc_qc_param = params_sc_min_quality(),
  streaming = TRUE,
  cell_id_col = NULL,
  batch_size = 1000L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  assertScMinQC(sc_qc_param)
  checkmate::qassert(streaming, "B1")
  checkmate::qassert(batch_size, "I1")
  checkmate::qassert(.verbose, "B1")

  # make rust happier
  h5_path <- path.expand(h5_path)

  # rust part
  h5_meta <- get_h5ad_dimensions(f_path = h5_path)

  rust_con <- get_sc_rust_ptr(object)

  file_res <- rust_con$h5_to_file(
    cs_type = h5_meta$type,
    h5_path = path.expand(h5_path),
    no_cells = h5_meta$dims["obs"],
    no_genes = h5_meta$dims["var"],
    qc_params = sc_qc_param,
    verbose = .verbose
  )

  if (streaming) {
    rust_con$generate_gene_based_data_streaming(
      batch_size = batch_size,
      verbose = .verbose
    )
  } else {
    rust_con$generate_gene_based_data(
      verbose = .verbose
    )
  }

  gene_nnz <- rust_con$get_nnz_genes(gene_indices = NULL)
  gene_nnz_dt <- data.table::data.table(no_cells_exp = gene_nnz)

  # duck db part
  duckdb_con <- get_sc_duckdb(object)
  if (.verbose) {
    message("Loading observations data from h5ad into the DuckDB.")
  }
  duckdb_con$populate_obs_from_h5(
    h5_path = h5_path,
    filter = as.integer(file_res$cell_indices + 1),
    cell_id_col = cell_id_col
  )
  if (.verbose) {
    message("Loading variables data from h5ad into the DuckDB.")
  }
  duckdb_con$populate_vars_from_h5(
    h5_path = h5_path,
    filter = as.integer(file_res$gene_indices + 1)
  )

  cell_res_dt <- data.table::setDT(file_res[c("nnz", "lib_size")])

  duckdb_con$add_data_obs(new_data = cell_res_dt)
  duckdb_con$add_data_var(new_data = gene_nnz_dt)
  duckdb_con$set_to_keep_column()
  cell_map <- duckdb_con$get_obs_index_map()
  gene_map <- duckdb_con$get_var_index_map()

  S7::prop(object, "dims") <- as.integer(rust_con$get_shape())
  object <- set_cell_mapping(x = object, cell_map = cell_map)
  object <- set_gene_mapping(x = object, gene_map = gene_map)

  return(object)
}

#### normalised ----------------------------------------------------------------

#' Load in h5ad with normalised counts to `SingleCells`
#'
#' @description
#' This function takes an h5ad file where only normalised counts are available
#' in the X slot and loads the obs and var data into the DuckDB of the
#' `SingleCells` class and the counts into a Rust-binarised format for
#' rapid access. Raw counts are reconstructed from the normalised values using
#' the library sizes stored in a specified obs column.
#'
#' The reconstruction assumes the normalisation was:
#' `norm = log1p(x / lib_size * target_size)`
#'
#' @param object `SingleCells` class.
#' @param h5_path File path to the h5ad object you wish to load in.
#' @param obs_lib_size_col String. Name of the obs column containing the total
#' counts per cell or spot (e.g. `"nCount_RNA"`). You will have to check the
#' original obs data for this.
#' @param target_size Numeric. The target size used in the original
#' normalisation (e.g. `1e4`).
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
#' @param cell_id_col Optional string. If a specific column in the h5ad obs
#' data represents the cell identifiers, you can specify it here.
#' @param streaming Boolean. Shall the data be streamed during the conversion
#' of CSR to CSC. Defaults to `TRUE` and should be used for larger data sets.
#' @param batch_size Integer. If `streaming = TRUE`, how many cells to process
#' in one batch. Defaults to `1000L`.
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @return It will populate the files on disk and return the class with updated
#' shape information.
#'
#' @export
load_h5ad_norm <- S7::new_generic(
  name = "load_h5ad_norm",
  dispatch_args = "object",
  fun = function(
    object,
    h5_path,
    obs_lib_size_col,
    target_size,
    sc_qc_param = params_sc_min_quality(),
    streaming = TRUE,
    cell_id_col = NULL,
    batch_size = 1000L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method load_h5ad_norm SingleCells
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(load_h5ad_norm, SingleCells) <- function(
  object,
  h5_path,
  obs_lib_size_col,
  target_size,
  sc_qc_param = params_sc_min_quality(),
  streaming = TRUE,
  cell_id_col = NULL,
  batch_size = 1000L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  assertScMinQC(sc_qc_param)
  checkmate::qassert(obs_lib_size_col, "S1")
  checkmate::qassert(target_size, "N1")
  checkmate::qassert(streaming, "B1")
  checkmate::qassert(batch_size, "I1")
  checkmate::qassert(.verbose, "B1")

  h5_path <- path.expand(h5_path)

  h5_meta <- get_h5ad_dimensions(f_path = h5_path)

  rust_con <- get_sc_rust_ptr(object)

  file_res <- rust_con$norm_h5_to_file(
    cs_type = h5_meta$type,
    h5_path = h5_path,
    no_cells = h5_meta$dims["obs"],
    no_genes = h5_meta$dims["var"],
    obs_lib_size_col = obs_lib_size_col,
    target_size = target_size,
    qc_params = sc_qc_param,
    verbose = .verbose
  )

  if (streaming) {
    rust_con$generate_gene_based_data_streaming(
      batch_size = batch_size,
      verbose = .verbose
    )
  } else {
    rust_con$generate_gene_based_data(
      verbose = .verbose
    )
  }

  gene_nnz <- rust_con$get_nnz_genes(gene_indices = NULL)
  gene_nnz_dt <- data.table::data.table(no_cells_exp = gene_nnz)

  duckdb_con <- get_sc_duckdb(object)
  if (.verbose) {
    message("Loading observations data from h5ad into the DuckDB.")
  }
  duckdb_con$populate_obs_from_h5(
    h5_path = h5_path,
    filter = as.integer(file_res$cell_indices + 1),
    cell_id_col = cell_id_col
  )
  if (.verbose) {
    message("Loading variables data from h5ad into the DuckDB.")
  }
  duckdb_con$populate_vars_from_h5(
    h5_path = h5_path,
    filter = as.integer(file_res$gene_indices + 1)
  )

  cell_res_dt <- data.table::setDT(file_res[c("nnz", "lib_size")])

  duckdb_con$add_data_obs(new_data = cell_res_dt)
  duckdb_con$add_data_var(new_data = gene_nnz_dt)
  duckdb_con$set_to_keep_column()
  cell_map <- duckdb_con$get_obs_index_map()
  gene_map <- duckdb_con$get_var_index_map()

  S7::prop(object, "dims") <- as.integer(rust_con$get_shape())
  object <- set_cell_mapping(x = object, cell_map = cell_map)
  object <- set_gene_mapping(x = object, gene_map = gene_map)

  return(object)
}

#### slower streaming version --------------------------------------------------

#' Stream in h5ad to `SingleCells`
#'
#' @description
#' This function takes an h5ad file and loads (via streaming) the obs and var
#' data into the DuckDB of the `SingleCells` class and the counts into
#' a Rust-binarised format for rapid access. During the reading in of the
#' counts, the log CPM transformation will occur automatically. This function
#' is specifically designed to deal with larger amounts of data and is slower
#' than [bixverse::load_h5ad()].
#'
#' @param object `SingleCells` class.
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
#' @param max_genes_in_memory Integer. How many genes shall be held in memory
#' at a given point. Defaults to `2000L`.
#' @param cell_batch_size Integer. How big are the batch sizes for the cells
#' in the transformation from the cell-based to gene-based format. Defaults to
#' `100000L`.
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @return It will populate the files on disk and return the class with updated
#' shape information.
#'
#' @export
stream_h5ad <- S7::new_generic(
  name = "stream_h5ad",
  dispatch_args = "object",
  fun = function(
    object,
    h5_path,
    sc_qc_param = params_sc_min_quality(),
    max_genes_in_memory = 2000L,
    cell_batch_size = 100000L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method stream_h5ad SingleCells
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(stream_h5ad, SingleCells) <- function(
  object,
  h5_path,
  sc_qc_param = params_sc_min_quality(),
  max_genes_in_memory = 2000L,
  cell_batch_size = 100000L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  assertScMinQC(sc_qc_param)
  checkmate::qassert(max_genes_in_memory, "I1")
  checkmate::qassert(max_genes_in_memory, "I1")
  checkmate::qassert(.verbose, "B1")

  # rust part
  h5_meta <- get_h5ad_dimensions(f_path = h5_path)

  rust_con <- get_sc_rust_ptr(object)

  file_res <- rust_con$h5_to_file_streaming(
    cs_type = h5_meta$type,
    h5_path = path.expand(h5_path),
    no_cells = h5_meta$dims["obs"],
    no_genes = h5_meta$dims["var"],
    qc_params = sc_qc_param,
    verbose = .verbose
  )

  rust_con$generate_gene_based_data_memory_bounded(
    max_genes_in_memory = max_genes_in_memory,
    cell_batch_size = cell_batch_size,
    verbose = .verbose
  )

  gene_nnz <- rust_con$get_nnz_genes(gene_indices = NULL)
  gene_nnz_dt <- data.table::data.table(no_cells_exp = gene_nnz)

  # duck db part
  duckdb_con <- get_sc_duckdb(object)
  if (.verbose) {
    message("Loading observations data from h5ad into the DuckDB.")
  }
  duckdb_con$populate_obs_from_h5(
    h5_path = h5_path,
    filter = as.integer(file_res$cell_indices + 1)
  )
  if (.verbose) {
    message("Loading variables data from h5ad into the DuckDB.")
  }
  duckdb_con$populate_vars_from_h5(
    h5_path = h5_path,
    filter = as.integer(file_res$gene_indices + 1)
  )

  cell_res_dt <- data.table::setDT(file_res[c("nnz", "lib_size")])

  duckdb_con$add_data_obs(new_data = cell_res_dt)
  duckdb_con$add_data_var(new_data = gene_nnz_dt)
  duckdb_con$set_to_keep_column()
  cell_map <- duckdb_con$get_obs_index_map()
  gene_map <- duckdb_con$get_var_index_map()

  S7::prop(object, "dims") <- as.integer(rust_con$get_shape())
  object <- set_cell_mapping(x = object, cell_map = cell_map)
  object <- set_gene_mapping(x = object, gene_map = gene_map)

  return(object)
}

##### multiple h5ad files ------------------------------------------------------

#' Load multiple h5ad files into a single `SingleCells`
#'
#' @description
#' Takes a pre-scan result from [bixverse::prescan_h5ad_files()] and loads
#' all files into a single experiment with global gene QC and sequential
#' cell indexing.
#'
#' @param object `SingleCells` class.
#' @param prescan_result Output of [bixverse::prescan_h5ad_files()].
#' @param sc_qc_param List. Output of [bixverse::params_sc_min_quality()].
#' @param cell_id_col Optional string. Column name for cell identifiers in obs.
#' @param streaming Boolean. Use streaming for CSR-to-CSC conversion.
#' @param batch_size Integer. Batch size if `streaming = TRUE`.
#' @param .verbose Boolean.
#'
#' @return The class with updated shape and populated DuckDB.
#'
#' @export
load_multi_h5ad <- S7::new_generic(
  name = "load_multi_h5ad",
  dispatch_args = "object",
  fun = function(
    object,
    prescan_result,
    sc_qc_param = params_sc_min_quality(),
    cell_id_col = NULL,
    streaming = TRUE,
    batch_size = 1000L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method load_multi_h5ad SingleCells
#'
#' @export
S7::method(load_multi_h5ad, SingleCells) <- function(
  object,
  prescan_result,
  sc_qc_param = params_sc_min_quality(),
  cell_id_col = NULL,
  streaming = TRUE,
  batch_size = 1000L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  assertScMinQC(sc_qc_param)
  checkmate::assertList(prescan_result)
  checkmate::assertTRUE(all(
    c("universe", "universe_size", "file_tasks") %in%
      names(prescan_result)
  ))
  checkmate::qassert(streaming, "B1")
  checkmate::qassert(batch_size, "I1")
  checkmate::qassert(.verbose, "B1")

  rust_con <- get_sc_rust_ptr(object)

  # call Rust
  file_res <- rust_con$multi_h5_to_file(
    file_tasks = prescan_result$file_tasks,
    universe_size = as.integer(prescan_result$universe_size),
    qc_params = sc_qc_param,
    verbose = .verbose
  )

  # generate gene-based binary
  if (streaming) {
    rust_con$generate_gene_based_data_streaming(
      batch_size = batch_size,
      verbose = .verbose
    )
  } else {
    rust_con$generate_gene_based_data(verbose = .verbose)
  }

  gene_nnz <- rust_con$get_nnz_genes(gene_indices = NULL)
  gene_nnz_dt <- data.table::data.table(no_cells_exp = gene_nnz)

  # DuckDB: obs
  duckdb_con <- get_sc_duckdb(object)

  if (.verbose) {
    message("Loading observation data from h5ad files into DuckDB.")
  }

  per_file_info <- lapply(file_res$per_file, function(f) {
    list(
      h5_path = prescan_result$file_tasks[[f$exp_id]]$h5_path,
      exp_id = f$exp_id,
      cell_filter = as.integer(f$cell_indices + 1L)
    )
  })

  duckdb_con$populate_obs_from_multi_h5(
    per_file_info = per_file_info,
    cell_id_col = cell_id_col
  )

  # DuckDB: var (use first file as reference)
  if (.verbose) {
    message("Loading variable data into DuckDB.")
  }

  # global_gene_indices maps final -> universe position (0-indexed)
  # universe genes are sorted alphabetically, so we can index directly
  global_gene_filter <- as.integer(file_res$global_gene_indices + 1L)

  # use first file for var metadata
  ref_h5_path <- prescan_result$file_tasks[[1L]]$h5_path
  ref_gene_filter <- which(
    !is.na(prescan_result$file_tasks[[1L]]$gene_local_to_universe)
  )
  # reorder to universe order, then filter to final set
  ref_universe_order <- order(
    prescan_result$file_tasks[[1L]]$gene_local_to_universe[ref_gene_filter]
  )
  ref_gene_filter <- ref_gene_filter[ref_universe_order][global_gene_filter]

  duckdb_con$populate_vars_from_h5(
    h5_path = ref_h5_path,
    filter = ref_gene_filter
  )

  # add QC columns
  per_file_qc <- lapply(file_res$per_file, function(f) {
    data.table::data.table(nnz = f$nnz, lib_size = f$lib_size)
  })
  cell_res_dt <- data.table::rbindlist(per_file_qc)

  duckdb_con$add_data_obs(new_data = cell_res_dt)
  duckdb_con$add_data_var(new_data = gene_nnz_dt)
  duckdb_con$set_to_keep_column()

  cell_map <- duckdb_con$get_obs_index_map()
  gene_map <- duckdb_con$get_var_index_map()

  S7::prop(object, "dims") <- as.integer(rust_con$get_shape())
  object <- set_cell_mapping(x = object, cell_map = cell_map)
  object <- set_gene_mapping(x = object, gene_map = gene_map)

  return(object)
}

### mtx ------------------------------------------------------------------------

#' Load in mtx/plain text files to `SingleCells`
#'
#' @description
#' This is a helper function to load in mtx files and corresponding plain text
#' files. It will automatically filter out low quality cells and only keep
#' high quality cells. Under the hood DucKDB and high performance Rust binary
#' files are being used to store the counts.
#'
#' @param object `SingleCells` class.
#' @param sc_mtx_io_param List. Please generate this one via
#' [bixverse::params_sc_mtx_io()]. Needs to contain:
#' \itemize{
#'   \item path_mtx - String. Path to the .mtx file
#'   \item path_obs - String. Path to the file containing cell/barcode info.
#'   \item path_var - String. String. Path to the file containing gene/variable
#'   info.
#'   \item cells_as_rows - Boolean. Do cells represent the rows or columns.
#' }
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
#' @param streaming Boolean. Shall the data be streamed during the conversion
#' of CSR to CSC. Defaults to `TRUE` and should be used for larger data sets.
#' @param batch_size Integer. If `streaming = TRUE`, how many cells to process
#' in one batch. Defaults to `1000L`.
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
    sc_mtx_io_param = params_sc_mtx_io(),
    sc_qc_param = params_sc_min_quality(),
    streaming = TRUE,
    batch_size = 1000L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method load_mtx SingleCells
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(load_mtx, SingleCells) <- function(
  object,
  sc_mtx_io_param = params_sc_mtx_io(),
  sc_qc_param = params_sc_min_quality(),
  streaming = TRUE,
  batch_size = 1000L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertClass(object, "bixverse::SingleCells")
  assertScMtxIO(sc_mtx_io_param)
  assertScMinQC(sc_qc_param)
  checkmate::qassert(.verbose, "B1")

  # rust part
  rust_con <- get_sc_rust_ptr(object)

  file_res <- with(
    sc_mtx_io_param,
    rust_con$mtx_to_file(
      mtx_path = path_mtx,
      qc_params = sc_qc_param,
      cells_as_rows = cells_as_rows,
      verbose = .verbose
    )
  )

  if (streaming) {
    rust_con$generate_gene_based_data_streaming(
      batch_size = batch_size,
      verbose = .verbose
    )
  } else {
    rust_con$generate_gene_based_data(
      verbose = .verbose
    )
  }

  gene_nnz <- rust_con$get_nnz_genes(gene_indices = NULL)
  gene_nnz_dt <- data.table::data.table(no_cells_exp = gene_nnz)

  # duckDB part
  duckdb_con <- get_sc_duckdb(object)

  with(
    sc_mtx_io_param,
    {
      if (.verbose) {
        message("Loading observations data from flat file into the DuckDB.")
      }
      duckdb_con$populate_obs_from_plain_text(
        f_path = path_obs,
        has_hdr = has_hdr,
        filter = as.integer(file_res$cell_indices + 1)
      )
      if (.verbose) {
        message("Loading variable data from flat file into the DuckDB.")
      }
      duckdb_con$populate_var_from_plain_text(
        f_path = path_var,
        has_hdr = has_hdr,
        filter = as.integer(file_res$gene_indices + 1)
      )
    }
  )

  cell_res_dt <- data.table::setDT(file_res[c("nnz", "lib_size")])

  duckdb_con$add_data_obs(new_data = cell_res_dt)
  duckdb_con$add_data_var(new_data = gene_nnz_dt)
  duckdb_con$set_to_keep_column()
  cell_map <- duckdb_con$get_obs_index_map()
  gene_map <- duckdb_con$get_var_index_map()

  S7::prop(object, "dims") <- as.integer(rust_con$get_shape())
  object <- set_cell_mapping(x = object, cell_map = cell_map)
  object <- set_gene_mapping(x = object, gene_map = gene_map)

  return(object)
}

### save to disk ---------------------------------------------------------------

#' Save memory-bound data to disk
#'
#' @description
#' Helper function that stores the memory-bound data to disk for checkpointing
#' or when you close the session for quick recovery of prior work. You have the
#' option to save as `".rds"` or `".qs2"` (you need to have the package `"qs2"`
#' installed for this option!).
#'
#' @param object `SingleCells` class.
#' @param type String. One of `c("qs2", "rds")`. Defines which binary format to
#' use. Will default to `"qs2"` for speed.
#'
#' @returns The object with added information on the data on disk.
#'
#' @export
save_sc_exp_to_disk <- S7::new_generic(
  name = "save_sc_exp_to_disk",
  dispatch_args = "object",
  fun = function(
    object,
    type = c("qs2", "rds")
  ) {
    S7::S7_dispatch()
  }
)

#' @method save_sc_exp_to_disk SingleCells
#'
#' @export
S7::method(save_sc_exp_to_disk, SingleCells) <- function(
  object,
  type = c("qs2", "rds")
) {
  type <- match.arg(type)
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::assertChoice(type, c("qs2", "rds"))

  # pull the data out from the class
  sc_map <- get_sc_map(object)
  sc_cache <- get_sc_cache(object)
  dir <- S7::prop(object, "dir_data")

  to_save <- list(sc_map = sc_map, sc_cache = sc_cache)

  if (type == "qs2") {
    if (!requireNamespace("qs2", quietly = TRUE)) {
      stop("Package 'qs2' is required to use qs2 format. Please install it.")
    }
    qs2::qs_save(to_save, file = file.path(dir, "memory.qs2"))
  } else if (type == "rds") {
    saveRDS(to_save, file = file.path(dir, "memory.rds"))
  }
}

### from disk ------------------------------------------------------------------

#' Load an existing SingleCells from disk
#'
#' @description
#' Helper function that can load the parameters to access the on-disk stored
#' data into the class.
#'
#' @param object `SingleCells` class.
#'
#' @returns The object with added information on the data on disk.
#'
#' @export
load_existing <- S7::new_generic(
  name = "load_existing",
  dispatch_args = "object",
  fun = function(
    object
  ) {
    S7::S7_dispatch()
  }
)

#' @method load_mtx SingleCells
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(load_existing, SingleCells) <- function(object) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))

  dir_data <- S7::prop(object, "dir_data")

  checkmate::assertFileExists(file.path(dir_data, "counts_cells.bin"))
  checkmate::assertFileExists(file.path(dir_data, "counts_genes.bin"))
  checkmate::assertFileExists(file.path(dir_data, "sc_duckdb.db"))

  # function body
  rust_con <- get_sc_rust_ptr(object)
  rust_con$set_from_file()

  duckdb_con <- get_sc_duckdb(object)

  if (any(c("memory.qs2", "memory.rds") %in% list.files(dir_data))) {
    message(paste(
      "Found stored data from save_sc_exp_to_disk().",
      "Loading that one into the object."
    ))

    # preferentially load qs2
    saved_data <- if ("memory.qs2" %in% list.files(dir_data)) {
      if (!requireNamespace("qs2", quietly = TRUE)) {
        stop("Package 'qs2' is required to use qs2 format. Please install it.")
      }
      qs2::qs_read(file.path(dir_data, "memory.qs2"))
    } else {
      readRDS(file.path(dir_data, "memory.rds"))
    }

    S7::prop(object, "dims") <- as.integer(rust_con$get_shape())
    S7::prop(object, "sc_map") <- saved_data$sc_map
    S7::prop(object, "sc_cache") <- saved_data$sc_cache

    # check that memory-stored data agrees with Rust to avoid panics...

    if (
      (length(get_cell_names(object)) != rust_con$get_shape()[1]) |
        (length(get_gene_names(object)) != rust_con$get_shape()[2])
    ) {
      stop(paste(
        "The dimensions of the found data do not agree with the Rust data."
      ))
    }
  } else {
    cell_map <- duckdb_con$get_obs_index_map()
    gene_map <- duckdb_con$get_var_index_map()
    cells_to_keep <- duckdb_con$get_cells_to_keep()
    S7::prop(object, "dims") <- as.integer(rust_con$get_shape())

    if (
      (length(cell_map) != rust_con$get_shape()[1]) |
        (length(gene_map) != rust_con$get_shape()[2])
    ) {
      stop(paste(
        "The data in the observation table and or var table do not match",
        "with what is stored on disk. Loading of the file failed"
      ))
    }

    object <- set_cell_mapping(x = object, cell_map = cell_map)
    object <- set_gene_mapping(x = object, gene_map = gene_map)
    S7::prop(object, "sc_map") <- set_cells_to_keep(
      S7::prop(object, "sc_map"),
      cells_to_keep
    )
  }

  return(object)
}
