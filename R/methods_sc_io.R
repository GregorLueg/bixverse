# single cell i/o --------------------------------------------------------------

## helpers ---------------------------------------------------------------------

### mtx ------------------------------------------------------------------------

#' Prescan multiple mtx directories for a multi-load
#'
#' @description
#' Walks each input directory, reads the features file to build the
#' **intersection** of gene IDs across inputs (matched by the first column,
#' typically Ensembl gene IDs), decompresses any `.mtx.gz` files into a
#' temporary directory, and builds the file tasks expected by
#' [bixverse::load_multi_mtx()].
#'
#' Each input directory must contain the standard 10x trio: a `.mtx` (or
#' `.mtx.gz`) file, a barcodes file, and a features/genes file. File names are
#' matched by extension; if a directory contains multiple matching files an
#' error is raised.
#'
#' @param dirs Character vector of input directories. Length >= 2.
#' @param exp_ids Character vector of experiment identifiers, one per directory.
#' Must be unique.
#' @param cells_as_rows Boolean. Applied uniformly to all inputs. Defaults
#' to `FALSE` (10x convention: genes are rows).
#' @param has_hdr Boolean. Whether the barcodes/features files have a
#' header row. Applied uniformly. Defaults to `FALSE` (10x convention).
#'
#' @return A list with:
#' \itemize{
#'   \item universe - Character vector of gene IDs in the intersection, in the
#'   order they will appear in the final var table.
#'   \item universe_size - Length of the universe.
#'   \item file_tasks - List of per-input task lists for Rust and DuckDB.
#'   \item temp_files - Character vector of temp files created during
#'   decompression; the caller should `unlink()` these after use.
#' }
#'
#' @export
prescan_mtx_dirs <- function(
  dirs,
  exp_ids,
  cells_as_rows = FALSE,
  has_hdr = FALSE
) {
  checkmate::assertCharacter(dirs, min.len = 2L)
  for (d in dirs) {
    checkmate::assertDirectoryExists(d)
  }
  checkmate::assertCharacter(exp_ids, len = length(dirs), unique = TRUE)
  checkmate::qassert(cells_as_rows, "B1")
  checkmate::qassert(has_hdr, "B1")

  temp_dir <- tempfile(pattern = "bixverse_mtx_prescan_")
  dir.create(temp_dir)
  temp_files <- character()

  # locate the trio inside each directory
  locate <- function(dir, pat) {
    files <- list.files(
      dir,
      pattern = pat,
      full.names = TRUE,
      ignore.case = TRUE
    )
    if (length(files) == 0L) {
      stop(sprintf("No file matching '%s' in %s", pat, dir))
    }
    if (length(files) > 1L) {
      stop(sprintf("Multiple files matching '%s' in %s", pat, dir))
    }
    files
  }

  # gunzip helper - returns decompressed path, registers temp file
  gunzip_to_temp <- function(path) {
    out_name <- sub("\\.gz$", "", basename(path), ignore.case = TRUE)
    out_path <- file.path(
      temp_dir,
      paste0(
        tools::file_path_sans_ext(out_name),
        "_",
        basename(tempfile("")),
        ".",
        tools::file_ext(out_name)
      )
    )
    con_in <- gzfile(path, open = "rb")
    on.exit(close(con_in), add = TRUE)
    con_out <- file(out_path, open = "wb")
    on.exit(close(con_out), add = TRUE)
    repeat {
      chunk <- readBin(con_in, "raw", n = 8 * 1024 * 1024)
      if (length(chunk) == 0L) {
        break
      }
      writeBin(chunk, con_out)
    }
    out_path
  }

  gene_sets <- vector("list", length(dirs))
  file_tasks <- vector("list", length(dirs))

  for (i in seq_along(dirs)) {
    d <- dirs[i]
    mtx_path <- locate(d, "\\.mtx(\\.gz)?$")
    features_path <- locate(d, "(features|genes)\\.(tsv|csv)(\\.gz)?$")
    barcodes_path <- locate(d, "barcodes\\.(tsv|csv)(\\.gz)?$")

    if (grepl("\\.gz$", mtx_path, ignore.case = TRUE)) {
      decompressed <- gunzip_to_temp(mtx_path)
      temp_files <- c(temp_files, decompressed)
      mtx_path <- decompressed
    }

    delim <- if (grepl("\\.tsv(\\.gz)?$", features_path, ignore.case = TRUE)) {
      "\t"
    } else {
      ","
    }
    feats <- data.table::fread(
      file = features_path,
      sep = delim,
      header = has_hdr
    )
    gene_ids <- as.character(feats[[1L]])

    gene_sets[[i]] <- gene_ids
    file_tasks[[i]] <- list(
      exp_id = exp_ids[i],
      mtx_path = mtx_path,
      barcodes_path = barcodes_path,
      features_path = features_path,
      cells_as_rows = cells_as_rows,
      has_hdr = has_hdr,
      local_gene_ids = gene_ids
    )
  }

  universe <- Reduce(intersect, gene_sets)
  if (length(universe) == 0L) {
    stop("Gene intersection across inputs is empty.")
  }

  for (i in seq_along(file_tasks)) {
    match_idx <- match(file_tasks[[i]]$local_gene_ids, universe)
    file_tasks[[i]]$gene_local_to_universe <- as.integer(match_idx - 1L)
    file_tasks[[i]]$local_gene_ids <- NULL
  }

  list(
    universe = universe,
    universe_size = length(universe),
    file_tasks = file_tasks,
    temp_files = temp_files
  )
}

### CSR to CSC conversion ------------------------------------------------------

#' Dispatch CSR-to-CSC generation based on streaming level
#'
#' Internal helper to keep the streaming dispatch consistent across all loaders.
#' Validates the streaming level and routes to the appropriate Rust method on
#' the count connector.
#'
#' @param rust_con The Rust count connector.
#' @param streaming Integer. `0L` for in-memory, `1L` for light streaming,
#' `2L` for heavy streaming with memory upper boundaries.
#' @param batch_size Integer. Batch size for light streaming.
#' @param max_genes_in_memory Integer. Genes held in memory at once for
#' heavy streaming.
#' @param cell_batch_size Integer. Cell batch size for heavy streaming.
#' @param .verbose Boolean.
#'
#' @return Invisible NULL. Side effect is the gene-based binary file.
#'
#' @keywords internal
.dispatch_gene_based_data <- function(
  rust_con,
  streaming,
  batch_size,
  max_genes_in_memory,
  cell_batch_size,
  .verbose
) {
  checkmate::qassert(streaming, "I1")
  checkmate::assertTRUE(streaming %in% c(0L, 1L, 2L))
  checkmate::qassert(batch_size, "I1")
  checkmate::qassert(max_genes_in_memory, "I1")
  checkmate::qassert(cell_batch_size, "I1")
  checkmate::qassert(.verbose, "B1")

  if (streaming == 0L) {
    if (.verbose) {
      message(" Loading data directly into memory for CSR to CSC conversion.")
    }
    rust_con$generate_gene_based_data(verbose = .verbose)
  } else if (streaming == 1L) {
    if (.verbose) {
      message(" Using light streaming for the CSR to CSC conversion.")
    }
    rust_con$generate_gene_based_data_streaming(
      batch_size = batch_size,
      verbose = .verbose
    )
  } else {
    if (.verbose) {
      message(paste(
        " Using heavy streaming with reduced memory",
        "pressure for the CSR to CSC conversion."
      ))
    }
    rust_con$generate_gene_based_data_memory_bounded(
      max_genes_in_memory = max_genes_in_memory,
      cell_batch_size = cell_batch_size,
      verbose = .verbose
    )
  }

  invisible(NULL)
}

## seurat ----------------------------------------------------------------------

#' Load in Seurat to `SingleCells`
#'
#' @description
#' This function takes a Seurat object and generates a `SingleCells` class
#' from it. The raw counts are extracted, written to the Rust binary format,
#' and the metadata is loaded into the DuckDB.
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
#' @param streaming Integer. CSR-to-CSC conversion mode. `0L` -> in-memory
#' (fastest, highest memory), `1L` -> light streaming with cell batching,
#' `2L` -> heavy streaming with memory upper boundaries. Defaults to `1L`.
#' @param batch_size Integer. Cell batch size when `streaming = 1L`. Defaults
#' to `1000L`.
#' @param max_genes_in_memory Integer. Maximum genes held in memory at once
#' when `streaming = 2L`. Defaults to `2000L`.
#' @param cell_batch_size Integer. Cell batch size when `streaming = 2L`.
#' Defaults to `100000L`.
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @return It will populate the files on disk and return the class with updated
#' shape information.
#'
#' @export
load_seurat <- S7::new_generic(
  name = "load_seurat",
  dispatch_args = "object",
  fun = function(
    object,
    seurat,
    sc_qc_param = params_sc_min_quality(),
    streaming = 1L,
    batch_size = 1000L,
    max_genes_in_memory = 2000L,
    cell_batch_size = 100000L,
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
  streaming = 1L,
  batch_size = 1000L,
  max_genes_in_memory = 2000L,
  cell_batch_size = 100000L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::assertClass(seurat, "Seurat")
  assertScMinQC(sc_qc_param)
  checkmate::qassert(streaming, "I1")
  checkmate::assertTRUE(streaming %in% c(0L, 1L, 2L))
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

  .dispatch_gene_based_data(
    rust_con = rust_con,
    streaming = streaming,
    batch_size = batch_size,
    max_genes_in_memory = max_genes_in_memory,
    cell_batch_size = cell_batch_size,
    .verbose = .verbose
  )

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
#' This function loads in data directly from R objects. The counts matrix must
#' be a `dgRMatrix` (rows = cells, columns = genes).
#'
#' @param object `SingleCells` class.
#' @param counts Sparse matrix. The cells represent the rows, the genes the
#' columns. Needs to be a `"dgRMatrix"`.
#' @param obs data.table. Cell metadata. Must have one row per cell in the
#' same order as `counts`.
#' @param var data.table. Feature metadata. Must have one row per gene in
#' the same order as `counts`.
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
#' @param streaming Integer. CSR-to-CSC conversion mode. `0L` -> in-memory
#' (fastest, highest memory), `1L` -> light streaming with cell batching,
#' `2L` -> heavy streaming with memory upper boundaries. Defaults to `1L`.
#' @param batch_size Integer. Cell batch size when `streaming = 1L`. Defaults
#' to `1000L`.
#' @param max_genes_in_memory Integer. Maximum genes held in memory at once
#' when `streaming = 2L`. Defaults to `2000L`.
#' @param cell_batch_size Integer. Cell batch size when `streaming = 2L`.
#' Defaults to `100000L`.
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @return It will populate the files on disk and return the class with updated
#' shape information.
#'
#' @export
load_r_data <- S7::new_generic(
  name = "load_r_data",
  dispatch_args = "object",
  fun = function(
    object,
    counts,
    obs,
    var,
    sc_qc_param = params_sc_min_quality(),
    streaming = 1L,
    batch_size = 1000L,
    max_genes_in_memory = 2000L,
    cell_batch_size = 100000L,
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
  streaming = 1L,
  batch_size = 1000L,
  max_genes_in_memory = 2000L,
  cell_batch_size = 100000L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::assertClass(counts, "dgRMatrix")
  no_cells <- nrow(counts)
  no_genes <- ncol(counts)
  checkmate::assertDataTable(obs, nrows = no_cells)
  checkmate::assertDataTable(var, nrows = no_genes)
  checkmate::qassert(streaming, "I1")
  checkmate::assertTRUE(streaming %in% c(0L, 1L, 2L))
  checkmate::qassert(.verbose, "B1")

  if (.verbose) {
    message("Writing counts to disk.")
  }

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

  .dispatch_gene_based_data(
    rust_con = rust_con,
    streaming = streaming,
    batch_size = batch_size,
    max_genes_in_memory = max_genes_in_memory,
    cell_batch_size = cell_batch_size,
    .verbose = .verbose
  )

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
#' @param cell_id_col Optional string. If a specific column in the h5ad obs
#' data is representing the cell identifiers, you can specify it here.
#' @param streaming Integer. `0L` -> all cells loaded in memory then transposed
#' (fastest, highest memory), `1L` -> light streaming with cell batching, `2L`
#' -> heavy streaming with memory upper boundaries on the gene side. Controls
#' memory pressure during the CSR-to-CSC conversion. Defaults to `1L`.
#' @param batch_size Integer. Cell batch size when `streaming = 1L`. Defaults
#' to `1000L`.
#' @param max_genes_in_memory Integer. Maximum genes held in memory at once
#' when `streaming = 2L`. Defaults to `2000L`.
#' @param cell_batch_size Integer. Cell batch size when `streaming = 2L`.
#' Defaults to `100000L`.
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
    streaming = 1L,
    cell_id_col = NULL,
    batch_size = 1000L,
    max_genes_in_memory = 2000L,
    cell_batch_size = 100000L,
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
  streaming = 1L,
  cell_id_col = NULL,
  batch_size = 1000L,
  max_genes_in_memory = 2000L,
  cell_batch_size = 100000L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  assertScMinQC(sc_qc_param)
  checkmate::qassert(streaming, "I1")
  checkmate::assertTRUE(streaming %in% c(0L, 1L, 2L))
  checkmate::qassert(.verbose, "B1")

  h5_path <- path.expand(h5_path)

  h5_meta <- get_h5ad_dimensions(f_path = h5_path)

  rust_con <- get_sc_rust_ptr(object)

  file_res <- rust_con$h5ad_to_file(
    cs_type = h5_meta$type,
    h5_path = h5_path,
    no_cells = h5_meta$dims["obs"],
    no_genes = h5_meta$dims["var"],
    qc_params = sc_qc_param,
    verbose = .verbose
  )

  .dispatch_gene_based_data(
    rust_con = rust_con,
    streaming = streaming,
    batch_size = batch_size,
    max_genes_in_memory = max_genes_in_memory,
    cell_batch_size = cell_batch_size,
    .verbose = .verbose
  )

  gene_nnz <- rust_con$get_nnz_genes(gene_indices = NULL)
  gene_nnz_dt <- data.table::data.table(no_cells_exp = gene_nnz)

  duckdb_con <- get_sc_duckdb(object)
  if (.verbose) {
    message("Loading observations data from h5ad into the DuckDB.")
  }
  duckdb_con$populate_obs_from_h5ad(
    h5_path = h5_path,
    filter = as.integer(file_res$cell_indices + 1),
    cell_id_col = cell_id_col
  )
  if (.verbose) {
    message("Loading variables data from h5ad into the DuckDB.")
  }
  duckdb_con$populate_vars_from_h5ad(
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
#' `SingleCells` class and the counts into a Rust-binarised format for rapid
#' access. Raw counts are reconstructed from the normalised values using the
#' library sizes stored in a specified obs column.
#'
#' The reconstruction assumes the normalisation was:
#' `norm = log1p(x / lib_size * target_size)`
#'
#' @param object `SingleCells` class.
#' @param h5_path File path to the h5ad object you wish to load in.
#' @param obs_lib_size_col String. Name of the obs column containing the total
#' counts per cell or spot (e.g. `"nCount_RNA"`).
#' @param target_size Numeric. The target size used in the original
#' normalisation (e.g. `1e4`).
#' @param sc_qc_param List. Output of [bixverse::params_sc_min_quality()].
#' @param cell_id_col Optional string. If a specific column in the h5ad obs
#' data represents the cell identifiers, you can specify it here.
#' @param streaming Integer. `0L` -> in-memory, `1L` -> light streaming, `2L`
#' -> heavy streaming with memory upper boundaries. Controls memory pressure
#' during CSR-to-CSC conversion. Defaults to `1L`.
#' @param batch_size Integer. Cell batch size when `streaming = 1L`. Defaults
#' to `1000L`.
#' @param max_genes_in_memory Integer. Maximum genes held in memory at once
#' when `streaming = 2L`. Defaults to `2000L`.
#' @param cell_batch_size Integer. Cell batch size when `streaming = 2L`.
#' Defaults to `100000L`.
#' @param .verbose Boolean.
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
    streaming = 1L,
    cell_id_col = NULL,
    batch_size = 1000L,
    max_genes_in_memory = 2000L,
    cell_batch_size = 100000L,
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
  streaming = 1L,
  cell_id_col = NULL,
  batch_size = 1000L,
  max_genes_in_memory = 2000L,
  cell_batch_size = 100000L,
  .verbose = TRUE
) {
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  assertScMinQC(sc_qc_param)
  checkmate::qassert(obs_lib_size_col, "S1")
  checkmate::qassert(target_size, "N1")
  checkmate::qassert(streaming, "I1")
  checkmate::assertTRUE(streaming %in% c(0L, 1L, 2L))
  checkmate::qassert(.verbose, "B1")

  h5_path <- path.expand(h5_path)

  h5_meta <- get_h5ad_dimensions(f_path = h5_path)

  rust_con <- get_sc_rust_ptr(object)

  file_res <- rust_con$norm_h5ad_to_file(
    cs_type = h5_meta$type,
    h5_path = h5_path,
    no_cells = h5_meta$dims["obs"],
    no_genes = h5_meta$dims["var"],
    obs_lib_size_col = obs_lib_size_col,
    target_size = target_size,
    qc_params = sc_qc_param,
    verbose = .verbose
  )

  .dispatch_gene_based_data(
    rust_con = rust_con,
    streaming = streaming,
    batch_size = batch_size,
    max_genes_in_memory = max_genes_in_memory,
    cell_batch_size = cell_batch_size,
    .verbose = .verbose
  )

  gene_nnz <- rust_con$get_nnz_genes(gene_indices = NULL)
  gene_nnz_dt <- data.table::data.table(no_cells_exp = gene_nnz)

  duckdb_con <- get_sc_duckdb(object)
  if (.verbose) {
    message("Loading observations data from h5ad into the DuckDB.")
  }
  duckdb_con$populate_obs_from_h5ad(
    h5_path = h5_path,
    filter = as.integer(file_res$cell_indices + 1),
    cell_id_col = cell_id_col
  )
  if (.verbose) {
    message("Loading variables data from h5ad into the DuckDB.")
  }
  duckdb_con$populate_vars_from_h5ad(
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

#' Stream in h5ad to `SingleCells` (alias)
#'
#' @description
#' Convenience alias for `load_h5ad(streaming = 2L)`. Kept for backwards
#' compatibility - forwards directly to [bixverse::load_h5ad()] with heavy
#' streaming enabled. Prefer calling `load_h5ad` directly with an explicit
#' `streaming` level.
#'
#' @param object `SingleCells` class.
#' @param h5_path File path to the h5ad object.
#' @param sc_qc_param List. Output of [bixverse::params_sc_min_quality()].
#' @param max_genes_in_memory Integer. Genes held in memory at once. Defaults
#' to `2000L`.
#' @param cell_batch_size Integer. Cell batch size. Defaults to `100000L`.
#' @param .verbose Boolean.
#'
#' @return The class with updated shape information.
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
S7::method(stream_h5ad, SingleCells) <- function(
  object,
  h5_path,
  sc_qc_param = params_sc_min_quality(),
  max_genes_in_memory = 2000L,
  cell_batch_size = 100000L,
  .verbose = TRUE
) {
  load_h5ad(
    object = object,
    h5_path = h5_path,
    sc_qc_param = sc_qc_param,
    streaming = 2L,
    max_genes_in_memory = max_genes_in_memory,
    cell_batch_size = cell_batch_size,
    .verbose = .verbose
  )
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
#' @param streaming Integer. `0L` -> in-memory, `1L` -> light streaming, `2L`
#' -> heavy streaming with memory upper boundaries. Defaults to `1L`.
#' @param batch_size Integer. Cell batch size when `streaming = 1L`. Defaults
#' to `1000L`.
#' @param max_genes_in_memory Integer. Maximum genes held in memory at once
#' when `streaming = 2L`. Defaults to `2000L`.
#' @param cell_batch_size Integer. Cell batch size when `streaming = 2L`.
#' Defaults to `100000L`.
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
    streaming = 1L,
    batch_size = 1000L,
    max_genes_in_memory = 2000L,
    cell_batch_size = 100000L,
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
  streaming = 1L,
  batch_size = 1000L,
  max_genes_in_memory = 2000L,
  cell_batch_size = 100000L,
  .verbose = TRUE
) {
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  assertScMinQC(sc_qc_param)
  checkmate::assertList(prescan_result)
  checkmate::assertTRUE(all(
    c("universe", "universe_size", "file_tasks") %in%
      names(prescan_result)
  ))
  checkmate::qassert(streaming, "I1")
  checkmate::assertTRUE(streaming %in% c(0L, 1L, 2L))
  checkmate::qassert(.verbose, "B1")

  rust_con <- get_sc_rust_ptr(object)

  file_res <- rust_con$multi_h5ad_to_file(
    file_tasks = prescan_result$file_tasks,
    universe_size = as.integer(prescan_result$universe_size),
    qc_params = sc_qc_param,
    verbose = .verbose
  )

  .dispatch_gene_based_data(
    rust_con = rust_con,
    streaming = streaming,
    batch_size = batch_size,
    max_genes_in_memory = max_genes_in_memory,
    cell_batch_size = cell_batch_size,
    .verbose = .verbose
  )

  gene_nnz <- rust_con$get_nnz_genes(gene_indices = NULL)
  gene_nnz_dt <- data.table::data.table(no_cells_exp = gene_nnz)

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

  duckdb_con$populate_obs_from_multi_h5ad(
    per_file_info = per_file_info,
    cell_id_col = cell_id_col
  )

  if (.verbose) {
    message("Loading variable data into DuckDB.")
  }

  final_gene_names <- prescan_result$universe[file_res$global_gene_indices + 1L]

  duckdb_con$populate_vars_from_h5ad_reordered(
    h5_path = prescan_result$file_tasks[[1L]]$h5_path,
    final_gene_names = final_gene_names
  )

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
#' [bixverse::params_sc_mtx_io()].
#' @param sc_qc_param List. Output of [bixverse::params_sc_min_quality()].
#' @param mtx_streaming Boolean. Shall the .mtx file ingestion itself be
#' streamed (via temp-file bucketing). Recommended for large mtx files.
#' Defaults to `TRUE`.
#' @param streaming Integer. CSR-to-CSC conversion mode. `0L` -> in-memory,
#' `1L` -> light streaming, `2L` -> heavy streaming with memory upper
#' boundaries. Defaults to `1L`.
#' @param batch_size Integer. Cell batch size when `streaming = 1L`. Defaults
#' to `1000L`.
#' @param max_genes_in_memory Integer. Maximum genes held in memory at once
#' when `streaming = 2L`. Defaults to `2000L`.
#' @param cell_batch_size Integer. Cell batch size when `streaming = 2L`.
#' Defaults to `100000L`.
#' @param .verbose Boolean.
#'
#' @return The class with updated shape information.
#'
#' @export
load_mtx <- S7::new_generic(
  name = "load_mtx",
  dispatch_args = "object",
  fun = function(
    object,
    sc_mtx_io_param = params_sc_mtx_io(),
    sc_qc_param = params_sc_min_quality(),
    mtx_streaming = TRUE,
    streaming = 1L,
    batch_size = 1000L,
    max_genes_in_memory = 2000L,
    cell_batch_size = 100000L,
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
  mtx_streaming = TRUE,
  streaming = 1L,
  batch_size = 1000L,
  max_genes_in_memory = 2000L,
  cell_batch_size = 100000L,
  .verbose = TRUE
) {
  checkmate::assertClass(object, "bixverse::SingleCells")
  assertScMtxIO(sc_mtx_io_param)
  assertScMinQC(sc_qc_param)
  checkmate::qassert(mtx_streaming, "B1")
  checkmate::qassert(streaming, "I1")
  checkmate::assertTRUE(streaming %in% c(0L, 1L, 2L))
  checkmate::qassert(.verbose, "B1")

  rust_con <- get_sc_rust_ptr(object)

  file_res <- if (mtx_streaming) {
    with(
      sc_mtx_io_param,
      rust_con$mtx_to_file_streaming(
        mtx_path = path_mtx,
        qc_params = sc_qc_param,
        cells_as_rows = cells_as_rows,
        verbose = .verbose
      )
    )
  } else {
    with(
      sc_mtx_io_param,
      rust_con$mtx_to_file(
        mtx_path = path_mtx,
        qc_params = sc_qc_param,
        cells_as_rows = cells_as_rows,
        verbose = .verbose
      )
    )
  }

  .dispatch_gene_based_data(
    rust_con = rust_con,
    streaming = streaming,
    batch_size = batch_size,
    max_genes_in_memory = max_genes_in_memory,
    cell_batch_size = cell_batch_size,
    .verbose = .verbose
  )

  gene_nnz <- rust_con$get_nnz_genes(gene_indices = NULL)
  gene_nnz_dt <- data.table::data.table(no_cells_exp = gene_nnz)

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

### multiple mtx ---------------------------------------------------------------

#' Load multiple mtx directories into a single `SingleCells`
#'
#' @description
#' Takes the result of [bixverse::prescan_mtx_dirs()] and loads all inputs
#' into a single experiment with global gene QC and sequential cell indexing.
#' The feature space is the **intersection** of input gene IDs.
#'
#' @param object `SingleCells` class.
#' @param prescan_result Output of [bixverse::prescan_mtx_dirs()].
#' @param sc_qc_param List. Output of [bixverse::params_sc_min_quality()].
#' @param streaming Integer. CSR-to-CSC conversion mode. `0L` -> in-memory,
#' `1L` -> light streaming, `2L` -> heavy streaming with memory upper
#' boundaries. Defaults to `1L`.
#' @param batch_size Integer. Cell batch size when `streaming = 1L`. Defaults
#' to `1000L`.
#' @param max_genes_in_memory Integer. Maximum genes held in memory at once
#' when `streaming = 2L`. Defaults to `2000L`.
#' @param cell_batch_size Integer. Cell batch size when `streaming = 2L`.
#' Defaults to `100000L`.
#' @param .verbose Boolean.
#'
#' @return The class with updated shape and populated DuckDB.
#'
#' @export
load_multi_mtx <- S7::new_generic(
  name = "load_multi_mtx",
  dispatch_args = "object",
  fun = function(
    object,
    prescan_result,
    sc_qc_param = params_sc_min_quality(),
    streaming = 1L,
    batch_size = 1000L,
    max_genes_in_memory = 2000L,
    cell_batch_size = 100000L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method load_multi_mtx SingleCells
#'
#' @export
S7::method(load_multi_mtx, SingleCells) <- function(
  object,
  prescan_result,
  sc_qc_param = params_sc_min_quality(),
  streaming = 1L,
  batch_size = 1000L,
  max_genes_in_memory = 2000L,
  cell_batch_size = 100000L,
  .verbose = TRUE
) {
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  assertScMinQC(sc_qc_param)
  checkmate::assertList(prescan_result)
  checkmate::assertTRUE(all(
    c("universe", "universe_size", "file_tasks", "temp_files") %in%
      names(prescan_result)
  ))
  checkmate::qassert(streaming, "I1")
  checkmate::assertTRUE(streaming %in% c(0L, 1L, 2L))
  checkmate::qassert(.verbose, "B1")

  on.exit(
    if (length(prescan_result$temp_files) > 0L) {
      unlink(prescan_result$temp_files)
    },
    add = TRUE
  )

  rust_con <- get_sc_rust_ptr(object)

  rust_tasks <- lapply(prescan_result$file_tasks, function(t) {
    list(
      exp_id = t$exp_id,
      mtx_path = t$mtx_path,
      cells_as_rows = t$cells_as_rows,
      gene_local_to_universe = t$gene_local_to_universe
    )
  })

  file_res <- rust_con$multi_mtx_to_file(
    file_tasks = rust_tasks,
    universe_size = as.integer(prescan_result$universe_size),
    qc_params = sc_qc_param,
    verbose = .verbose
  )

  .dispatch_gene_based_data(
    rust_con = rust_con,
    streaming = streaming,
    batch_size = batch_size,
    max_genes_in_memory = max_genes_in_memory,
    cell_batch_size = cell_batch_size,
    .verbose = .verbose
  )

  gene_nnz <- rust_con$get_nnz_genes(gene_indices = NULL)
  gene_nnz_dt <- data.table::data.table(no_cells_exp = gene_nnz)

  duckdb_con <- get_sc_duckdb(object)

  if (.verbose) {
    message("Loading barcodes from input directories into DuckDB.")
  }
  per_file_obs <- lapply(seq_along(prescan_result$file_tasks), function(i) {
    t <- prescan_result$file_tasks[[i]]
    list(
      f_path = t$barcodes_path,
      exp_id = t$exp_id,
      has_hdr = t$has_hdr,
      cell_filter = as.integer(file_res$per_file[[i]]$cell_indices + 1L)
    )
  })
  duckdb_con$populate_obs_from_multi_plain_text(per_file_info = per_file_obs)

  if (.verbose) {
    message("Loading features into DuckDB.")
  }
  final_gene_names <-
    prescan_result$universe[file_res$global_gene_indices + 1L]
  duckdb_con$populate_vars_from_plain_text_reordered(
    f_path = prescan_result$file_tasks[[1L]]$features_path,
    has_hdr = prescan_result$file_tasks[[1L]]$has_hdr,
    final_gene_names = final_gene_names
  )

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
