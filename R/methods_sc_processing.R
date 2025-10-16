# single cell processing methods -----------------------------------------------

## i/o -------------------------------------------------------------------------

### seurat ---------------------------------------------------------------------

#' Load in Seurat to `single_cell_exp`
#'
#' @description
#' This function takes a Seurat file and generates `single_cell_exp` class
#' from it.
#'
#' @param object `single_cell_exp` class.
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

#' @method load_seurat single_cell_exp
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(load_seurat, single_cell_exp) <- function(
  object,
  seurat,
  sc_qc_param = params_sc_min_quality(),
  batch_size = 1000L,
  streaming = TRUE,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, single_cell_exp))
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
  duckdb_con$set_to_keep_column()
  cell_map <- duckdb_con$get_obs_index_map()
  gene_map <- duckdb_con$get_var_index_map()

  S7::prop(object, "dims") <- as.integer(rust_con$get_shape())
  object <- set_cell_mapping(x = object, cell_map = cell_map)
  object <- set_gene_mapping(x = object, gene_map = gene_map)

  return(object)
}

### direct r -------------------------------------------------------------------

#' Load in Seurat to `single_cell_exp`
#'
#' @description
#' This function takes a Seurat file and generates `single_cell_exp` class
#' from it.
#'
#' @param object `single_cell_exp` class.
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

#' @method load_r_data single_cell_exp
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(load_r_data, single_cell_exp) <- function(
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
  checkmate::assertTRUE(S7::S7_inherits(object, single_cell_exp))
  checkmate::assertClass(counts, "dgRMatrix")
  no_cells <- nrow(counts)
  no_genes <- ncol(counts)
  checkmate::assertDataTable(obs, nrows = no_cells)
  checkmate::assertDataTable(var, nrows = no_genes)
  checkmate::qassert(.verbose, "B1")

  # body
  counts <- sparse_mat_to_list(counts)

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

  duckdb_con$add_data_obs(new_data = cell_res_dt)
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

#' Load in h5ad to `single_cell_exp`
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
    batch_size = 1000L,
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
  streaming = TRUE,
  batch_size = 1000L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, single_cell_exp))
  assertScMinQC(sc_qc_param)
  checkmate::qassert(streaming, "B1")
  checkmate::qassert(batch_size, "I1")
  checkmate::qassert(.verbose, "B1")

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
  duckdb_con$set_to_keep_column()
  cell_map <- duckdb_con$get_obs_index_map()
  gene_map <- duckdb_con$get_var_index_map()

  S7::prop(object, "dims") <- as.integer(rust_con$get_shape())
  object <- set_cell_mapping(x = object, cell_map = cell_map)
  object <- set_gene_mapping(x = object, gene_map = gene_map)

  return(object)
}

#### slower streaming version --------------------------------------------------

#' Stream in h5ad to `single_cell_exp`
#'
#' @description
#' This function takes an h5ad file and loads (via streaming) the obs and var
#' data into the DuckDB of the `single_cell_exp` class and the counts into
#' a Rust-binarised format for rapid access. During the reading in of the
#' counts, the log CPM transformation will occur automatically. This function
#' is specifically designed to deal with larger amounts of data and is slower
#' than [bixverse::load_h5ad()].
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

#' @method stream_h5ad single_cell_exp
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(stream_h5ad, single_cell_exp) <- function(
  object,
  h5_path,
  sc_qc_param = params_sc_min_quality(),
  max_genes_in_memory = 2000L,
  cell_batch_size = 100000L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, single_cell_exp))
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
  duckdb_con$set_to_keep_column()
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

#' @method load_mtx single_cell_exp
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(load_mtx, single_cell_exp) <- function(
  object,
  sc_mtx_io_param = params_sc_mtx_io(),
  sc_qc_param = params_sc_min_quality(),
  streaming = TRUE,
  batch_size = 1000L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertClass(object, "bixverse::single_cell_exp")
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
  duckdb_con$set_to_keep_column()
  cell_map <- duckdb_con$get_obs_index_map()
  gene_map <- duckdb_con$get_var_index_map()

  S7::prop(object, "dims") <- as.integer(rust_con$get_shape())
  object <- set_cell_mapping(x = object, cell_map = cell_map)
  object <- set_gene_mapping(x = object, gene_map = gene_map)

  return(object)
}

### from disk ------------------------------------------------------------------

#' Load an existing single_cell_exp from disk
#'
#' @description
#' Helper function that can load the parameters to access the on-disk stored
#' data into the class.
#'
#' @param object `single_cell_exp` class.
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

#' @method load_mtx single_cell_exp
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(load_existing, single_cell_exp) <- function(object) {
  dir_data <- S7::prop(object, "dir_data")

  checkmate::assertFileExists(file.path(dir_data, "counts_cells.bin"))
  checkmate::assertFileExists(file.path(dir_data, "counts_genes.bin"))
  checkmate::assertFileExists(file.path(dir_data, "sc_duckdb.db"))

  # function body
  rust_con <- get_sc_rust_ptr(object)
  rust_con$set_from_file()

  duckdb_con <- get_sc_duckdb(object)

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

  return(object)
}

## qc --------------------------------------------------------------------------

### doublet detection ----------------------------------------------------------

#' Doublet detection with Scrublet
#'
#' @description This function implements the doublet detection from Scrublet,
#' see Wolock, et al. Briefly, arteficial doublets are being generated from
#' the data via random combination of initial cells. Highly variable genes (HVG)
#' are being identified and the observed cells are being projected on a PCA
#' space. Subsequently, the simulated doublets are being projected on the same
#' PCA space given the same HVGs; kNN graphs are being generated and a kNN
#' classifier is used to assign a probability that a given cell in the original
#' data is a doublet. For more details, please check the publication.
#'
#' @param object `single_cell_exp` class.
#' @param target_size Numeric. The library target size for the simulated cells.
#' Needs to be the same as the initial one (usually `1e5`).
#' @param scrublet_params A list with the scrublet parameters, see
#' [bixverse::params_scrublet()]. This list contains:
#' \itemize{
#'   \item min_gene_var_pctl - Numeric. Percentile threshold for highly variable
#'   genes. For example, 0.85 means keep genes in top 15% of variability.
#'   \item hvg_method - String. Method for highly variable gene selection. One
#'   of c("vst", "mvb", "dispersion"). Defaults to "vst" (variance stabilising
#'   transformation).
#'   \item loess_span - Numeric. Span parameter for loess fitting in VST method.
#'   Controls smoothness of the fitted curve.
#'   \item clip_max - Optional numeric. Optional maximum value for clipping in
#'   variance stabilisation.
#'   \item sim_doublet_ratio - Numeric. Number of doublets to simulate relative
#'   to the number of observed cells. For example, 2.0 simulates twice as many
#'   doublets as there are cells. Defaults to `1.0`.
#'   \item expected_doublet_rate - Numeric. Expected doublet rate for the
#'   experiment, typically 0.05-0.10 depending on cell loading. Must be between
#'   0 and 1.
#'   \item stdev_doublet_rate - Numeric. Uncertainty in the expected doublet
#'   rate.
#'   \item no_pcs - Integer. Number of principal components to use for
#'   embedding.
#'   \item random_svd - Boolean. Whether to use randomised SVD (faster) vs
#'   exact SVD.
#'   \item k - Integer. Number of nearest neighbours for the kNN graph. If 0
#'   (default), automatically calculated as `round(0.5 * sqrt(n_cells))`.
#'   \item knn_method - String. Distance metric to use. One of
#'   `c("cosine", "euclidean")`. Defaults to `"cosine"`.
#'   \item search_budget - Integer. Search budget for Annoy algorithm (higher =
#'   more accurate but slower).
#'   \item n_trees - Integer. Number of trees for Annoy index generation.
#'   \item n_bins - Integer. Number of bins for histogram-based automatic
#'   threshold detection. Typically 50-100.
#'   \item manual_threshold - Optional numeric. Optional manual doublet score
#'   threshold. If NULL (default), threshold is automatically detected from
#'   simulated doublet score distribution.
#' }
#' @param seed Integer. Random seed.
#' @param streaming Boolean. Shall streaming be used during the HVG
#' calculations. Slower, but less memory usage.
#' @param return_combined_pca Boolean. Shall the PCA of the observed cells and
#' simulated doublets be returned.
#' @param return_pairs Boolean. Shall the pairs be returned.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @return A `scrublet_res` class that has with the following items:
#' \itemize{
#'   \item predicted_doublets - Boolean vector indicating which observed cells
#'   predicted as doublets (TRUE = doublet, FALSE = singlet).
#'   \item doublet_scores_obs - Numerical vector with the likelihood of being
#'   a doublet for the observed cells.
#'   \item doublet_scores_sim - Numerical vector with the likelihood of being
#'   a doublet for the simulated cells.
#'   \item doublet_errors_obs - Numerical vector with the standard errors of
#'   the scores for the observed cells.
#'   \item z_scores - Z-scores for the observed cells. Represents:
#'   `score - threshold / error`.
#'   \item threshold - Used threshold.
#'   \item detected_doublet_rate - Fraction of cells that are called as
#'   doublet.
#'   \item detectable_doublet_fraction - Fraction of simulated doublets with
#'   scores above the threshold.
#'   \item overall_doublet_rate - Estimated overall doublet rate. Should roughly
#'   match the expected doublet rate.
#' }
#'
#' @export
#'
#' @references Wollock, et al., Cell Syst, 2020
scrublet_sc <- S7::new_generic(
  name = "scrublet_sc",
  dispatch_args = "object",
  fun = function(
    object,
    scrublet_params = params_scrublet(),
    seed = 42L,
    streaming = FALSE,
    return_combined_pca = FALSE,
    return_pairs = FALSE,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method scrublet_sc single_cell_exp
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(scrublet_sc, single_cell_exp) <- function(
  object,
  scrublet_params = params_scrublet(),
  seed = 42L,
  streaming = FALSE,
  return_combined_pca = FALSE,
  return_pairs = FALSE,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, single_cell_exp))
  assertScScrublet(scrublet_params)
  checkmate::qassert(seed, "I1")
  checkmate::qassert(streaming, "B1")
  checkmate::qassert(return_combined_pca, "B1")
  checkmate::qassert(return_pairs, "B1")
  checkmate::qassert(.verbose, "B1")

  # function body
  cells_to_keep <- get_cells_to_keep(object)

  scrublet_res <- rs_sc_scrublet(
    f_path_gene = get_rust_count_gene_f_path(object),
    f_path_cell = get_rust_count_cell_f_path(object),
    cells_to_keep = cells_to_keep,
    scrublet_params = scrublet_params,
    seed = seed,
    verbose = .verbose,
    streaming = streaming,
    return_combined_pca = return_combined_pca,
    return_pairs = return_pairs
  )

  attr(scrublet_res, "cell_indices") <- cells_to_keep
  class(scrublet_res) <- "scrublet_res"

  return(scrublet_res)
}

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
#' @param streaming Boolean. Shall the cells be streamed in. Useful for larger
#' data sets where you wish to avoid loading in the whole data. Default to
#' `FALSE`.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @return It will add the columns based on the names in the `gene_set_list` to
#' the obs table.
#'
#' @export
gene_set_proportions_sc <- S7::new_generic(
  name = "gene_set_proportions_sc",
  dispatch_args = "object",
  fun = function(
    object,
    gene_set_list,
    streaming = FALSE,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method gene_set_proportions_sc single_cell_exp
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(gene_set_proportions_sc, single_cell_exp) <- function(
  object,
  gene_set_list,
  streaming = FALSE,
  .verbose = TRUE
) {
  checkmate::assertClass(object, "bixverse::single_cell_exp")
  checkmate::assertList(gene_set_list, names = "named", types = "character")
  checkmate::qassert(streaming, "B1")
  checkmate::qassert(.verbose, "B1")

  gene_set_list_tidy <- purrr::map(gene_set_list, \(g) {
    get_gene_indices(object, gene_ids = g, rust_index = TRUE)
  })
  names(gene_set_list_tidy) <- names(gene_set_list)

  rs_results <- rs_sc_get_gene_set_perc(
    f_path_cell = get_rust_count_cell_f_path(object),
    cell_indices = get_cells_to_keep(object),
    gene_set_idx = gene_set_list_tidy,
    streaming = streaming,
    verbose = .verbose
  )

  class(rs_results) <- "sc_proportion_res"
  attr(rs_results, "cell_indices") <- get_cells_to_keep(object)

  duckdb_con <- get_sc_duckdb(object)

  duckdb_con$join_data_obs(get_obs_data(rs_results))

  if (length(rs_results) == 1) {
    # need to deal with the case of only one gene set here...
    object[[names(rs_results)]] <- rs_results[[1]]
  } else {
    object[[names(rs_results)]] <- rs_results
  }

  return(object)
}

### hvg ------------------------------------------------------------------------

#' Identify HVGs
#'
#' @description
#' This is a helper function to identify highly variable genes. At the moment
#' the implementation has only the VST-based version (known as Seurat v3). The
#' other methods will be implemented in the future.
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
#' @param streaming Boolean. Shall the genes be streamed in. Useful for larger
#' data sets where you wish to avoid loading in the whole data. Defaults to
#' `FALSE`.
#' @param .verbose Boolean. Controls verbosity and returns run times.
#'
#' @return It will add the columns based on the names in the `gene_set_list` to
#' the obs table.
#'
#' @export
find_hvg_sc <- S7::new_generic(
  name = "find_hvg_sc",
  dispatch_args = "object",
  fun = function(
    object,
    hvg_no = 2000L,
    hvg_params = params_sc_hvg(),
    streaming = FALSE,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method find_hvg_sc single_cell_exp
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(find_hvg_sc, single_cell_exp) <- function(
  object,
  hvg_no = 2000L,
  hvg_params = params_sc_hvg(),
  streaming = FALSE,
  .verbose = TRUE
) {
  checkmate::assertClass(object, "bixverse::single_cell_exp")
  checkmate::qassert(hvg_no, "I1")
  assertScHvg(hvg_params)
  checkmate::qassert(streaming, "B1")
  checkmate::qassert(.verbose, "B1")

  if (length(get_cells_to_keep(object)) == 0) {
    warning(paste(
      "You need to set the cells to keep with set_cells_to_keep().",
      "Returning class as is."
    ))
    return(object)
  }

  res <- with(
    hvg_params,
    rs_sc_hvg(
      f_path_gene = get_rust_count_gene_f_path(object),
      hvg_method = method,
      cell_indices = get_cells_to_keep(object),
      loess_span = loess_span,
      clip_max = NULL,
      streaming = streaming,
      verbose = .verbose
    )
  )

  object <- set_sc_new_var_cols(object = object, data_list = res)

  hvg <- order(res$var_std, decreasing = TRUE)[1:hvg_no]

  object <- set_hvg(object, hvg = hvg)

  return(object)
}

## dimension reduction and knn/snn ---------------------------------------------

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
#' @return The function will add the PCA factors, loadings and singular values
#' to the object cache in memory.
#'
#' @export
calculate_pca_sc <- S7::new_generic(
  name = "calculate_pca_sc",
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

#' @method calculate_pca_sc single_cell_exp
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(calculate_pca_sc, single_cell_exp) <- function(
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
      "No HVGs identified in this object. Did you run find_hvg_sc()?",
      "Returning object as is."
    ))
    return(object)
  }

  zeallot::`%<-%`(
    c(pca_factors, pca_loadings, singular_values, scaled),
    rs_sc_pca(
      f_path_gene = get_rust_count_gene_f_path(object),
      no_pcs = no_pcs,
      random_svd = randomised_svd,
      cell_indices = get_cells_to_keep(object),
      gene_indices = get_hvg(object),
      seed = seed,
      return_scaled = FALSE,
      verbose = .verbose
    )
  )

  object <- set_pca_factors(object, pca_factors)
  object <- set_pca_loadings(object, pca_loadings)
  object <- set_pca_singular_vals(object, singular_values[1:no_pcs])

  return(object)
}

### neighbours -----------------------------------------------------------------

#' Find the neighbours for single cell.
#'
#' @description
#' This function will generate the kNNs based on a given embedding (atm,
#' only option is PCA). Two different algorithms are implemented with different
#' speed and accuracy to approximate the nearest neighbours. `"annoy"` is more
#' rapid and based on the `Approximate Nearest Neigbours Oh Yeah` algorithm,
#' whereas `"hnsw"` implements a `Hierarchical Navigatable Small Worlds` vector
#' search that is slower, but more precise. Subsequently, the kNN data will
#' be used to generate an sNN igraph for clustering methods.
#'
#' @param object `single_cell_exp` class.
#' @param embd_to_use String. The embedding to use. Atm, the only option is
#' `"pca"`.
#' @param no_embd_to_use Optional integer. Number of embedding dimensions to
#' use. If `NULL` all will be used.
#' @param neighbours_params List. Output of [bixverse::params_sc_neighbours()].
#' A list with the following items:
#' \itemize{
#'   \item k - Integer. Number of neighbours to identify.
#'   \item n_trees -  Integer. Number of trees to use for the `annoy` algorithm.
#'   The higher, the longer the algorithm takes, but the more precise the
#'   approximated nearest neighbours.
#'   \item search_budget - Integer. Search budget per tree for the `annoy`
#'   algorithm. The higher, the longer the algorithm takes, but the more precise
#'   the approximated nearest neighbours.
#'   \item knn_algorithm - String. One of `c("annoy", "hnsw")`. `"hnsw"` takes
#'   longer, is more precise and more memory friendly. `"annoy"` is faster, less
#'   precise and will take more memory.
#'   \item ann_dist - String. One of `c("cosine", "euclidean")`.
#'   \item full_snn - Boolean. Shall the sNN graph be generated across all
#'   cells (standard in the `bluster` package.) Defaults to `FALSE`.
#'   \item pruning - Value below which the weight in the sNN graph is set to 0.
#'   \item snn_similarity - String. One of `c("rank", "jaccard")`. Defines how
#'   the weight form the SNN graph is calculated. For details, please see
#'   [bixverse::params_sc_neighbours()].
#' }
#' @param seed Integer. For reproducibility.
#' @param .verbose Boolean. Controls verbosity and returns run times.
#'
#' @return The object with added KNN matrix.
#'
#' @export
find_neighbours_sc <- S7::new_generic(
  name = "find_neighbours_sc",
  dispatch_args = "object",
  fun = function(
    object,
    embd_to_use = "pca",
    no_embd_to_use = NULL,
    neighbours_params = params_sc_neighbours(),
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method find_neighbours_sc single_cell_exp
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(find_neighbours_sc, single_cell_exp) <- function(
  object,
  embd_to_use = "pca",
  no_embd_to_use = NULL,
  neighbours_params = params_sc_neighbours(),
  seed = 42L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertClass(object, "bixverse::single_cell_exp")
  checkmate::assertChoice(embd_to_use, c("pca"))
  checkmate::qassert(no_embd_to_use, c("I1", "0"))
  assertScNeighbours(neighbours_params)
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, "B1")

  # get the embedding
  embd <- switch(embd_to_use, pca = get_pca_factors(object))

  # early return
  if (is.null(embd)) {
    warning(
      paste(
        "The desired embedding was not found. Please check what you are doing",
        "Returning class as is"
      )
    )

    return(object)
  }

  if (!is.null(no_embd_to_use)) {
    to_take <- min(c(no_embd_to_use, ncol(embd)))
    embd <- embd[, 1:to_take]
  }

  if (.verbose) {
    message(sprintf(
      "Generating kNN data with %s method.",
      neighbours_params$knn_algorithm
    ))
  }

  knn_data <- with(
    neighbours_params,
    rs_sc_knn(
      embd = embd,
      no_neighbours = k,
      seed = seed,
      n_trees = n_trees,
      search_budget = search_budget,
      verbose = .verbose,
      algorithm_type = knn_algorithm,
      ann_dist = ann_dist
    )
  )

  object <- set_knn(object, knn_data)

  if (.verbose) {
    message(sprintf(
      "Generating sNN graph (full: %s).",
      neighbours_params$full_snn
    ))
  }

  snn_graph_rs <- with(
    neighbours_params,
    rs_sc_snn(
      knn_data,
      snn_method = snn_similarity,
      pruning = pruning,
      limited_graph = !full_snn,
      verbose = .verbose
    )
  )

  if (.verbose) {
    message("Transforming sNN data to igraph.")
  }

  snn_g <- igraph::make_graph(snn_graph_rs$edges + 1, directed = FALSE)
  igraph::E(snn_g)$weight <- snn_graph_rs$weights

  object <- set_snn_graph(object, snn_graph = snn_g)

  return(object)
}

### clustering -----------------------------------------------------------------

#' Graph-based clustering of cells on the sNN graph
#'
#' @description
#' This function will apply Leiden clustering on the sNN graph with the
#' given resolution and add a column to the obs table.
#'
#' @param object `single_cell_exp` class.
#' @param res Numeric. The resolution parameter for [igraph::cluster_leiden()].
#' @param name String. The name to add to the obs table in the DuckDB.
#'
#' @return The object with added clustering in the obs table.
#'
#' @export
find_clusters_sc <- S7::new_generic(
  name = "find_clusters_sc",
  dispatch_args = "object",
  fun = function(
    object,
    res = 1,
    name = "leiden_clustering"
  ) {
    S7::S7_dispatch()
  }
)

#' @method find_clusters_sc single_cell_exp
#'
#' @export
#'
#' @importFrom zeallot `%<-%`
#' @importFrom magrittr `%>%`
S7::method(find_clusters_sc, single_cell_exp) <- function(
  object,
  res = 1,
  name = "leiden_clustering"
) {
  # checks
  checkmate::assertClass(object, "bixverse::single_cell_exp")
  checkmate::qassert(res, "N1")
  checkmate::qassert(name, "S1")

  snn_graph <- get_snn_graph(object)

  if (is.null(snn_graph)) {
    warning(paste(
      "No sNN graph found. Did you run find_neighbours_sc()",
      "Returning class as is."
    ))
    return(object)
  }

  leiden_clusters <- igraph::cluster_leiden(
    snn_graph,
    objective_function = "modularity",
    resolution = res
  )

  duckdb_con <- get_sc_duckdb(object)

  new_data <- data.table::data.table(
    cell_idx = get_cells_to_keep(object) + 1, # needs to be 1-indexed
    new_data = leiden_clusters$membership
  )
  data.table::setnames(new_data, "new_data", name)

  duckdb_con$join_data_obs(new_data = new_data)

  return(object)
}
