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
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @returns A list with:
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
#'
#' @examples
#' # two CellRanger style directories reduced to their shared gene universe
#' data <- generate_single_cell_test_data(
#'   syn_data_params = params_sc_synthetic_data(n_cells = 200L, n_genes = 40L)
#' )
#' dirs <- c(tempfile("cr_a"), tempfile("cr_b"))
#' for (d in dirs) {
#'   dir.create(d, recursive = TRUE)
#'   write_cellranger_output(
#'     d, data$counts, data$obs, data$var,
#'     rows = "cells", format_type = "csv", .verbose = FALSE
#'   )
#' }
#' scan_res <- prescan_mtx_dirs(
#'   dirs = dirs,
#'   exp_ids = c("a", "b"),
#'   cells_as_rows = TRUE,
#'   has_hdr = TRUE,
#'   .verbose = FALSE
#' )
#' scan_res$universe_size
#'
#' unlink(c(dirs, scan_res$temp_files), recursive = TRUE, force = TRUE)
prescan_mtx_dirs <- function(
  dirs,
  exp_ids,
  cells_as_rows = FALSE,
  has_hdr = FALSE,
  .verbose = TRUE
) {
  checkmate::assertCharacter(dirs, min.len = 2L)
  for (d in dirs) {
    checkmate::assertDirectoryExists(d)
  }
  checkmate::assertCharacter(exp_ids, len = length(dirs), unique = TRUE)
  checkmate::qassert(cells_as_rows, "B1")
  checkmate::qassert(has_hdr, "B1")
  checkmate::qassert(.verbose, "B1")

  temp_dir <- tempfile(pattern = "bixverse_mtx_prescan_")
  dir.create(temp_dir)
  temp_files <- character()

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

  if (.verbose) {
    cli::cli_progress_bar("Prescanning mtx directories", total = length(dirs))
  }

  for (i in seq_along(dirs)) {
    if (.verbose) {
      cli::cli_progress_update(status = exp_ids[i])
    }

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

  if (.verbose) {
    cli::cli_progress_done()
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

#' Generate the gene-based binary from the cell-based one
#'
#' Internal helper so every loader runs the CSR to CSC conversion, and warns
#' about the retired streaming arguments, the same way.
#'
#' @param rust_con The Rust count connector.
#' @param csc_mem_gb Optional numeric. Memory in GB for the conversion buffers.
#' `NULL` converts in one pass.
#' @param streaming,batch_size,max_genes_in_memory,cell_batch_size Retired
#' arguments, forwarded from the caller only to warn if they were supplied.
#' @param .verbose Boolean.
#'
#' @returns Invisible NULL. Side effect is the gene-based binary file.
#'
#' @keywords internal
.dispatch_gene_based_data <- function(
  rust_con,
  csc_mem_gb,
  streaming = deprecated(),
  batch_size = deprecated(),
  max_genes_in_memory = deprecated(),
  cell_batch_size = deprecated(),
  .verbose
) {
  checkmate::qassert(csc_mem_gb, c("0", "N1(0,)"))
  checkmate::qassert(.verbose, "B1")

  supplied <- c(
    streaming = lifecycle::is_present(streaming),
    batch_size = lifecycle::is_present(batch_size),
    max_genes_in_memory = lifecycle::is_present(max_genes_in_memory),
    cell_batch_size = lifecycle::is_present(cell_batch_size)
  )
  if (any(supplied)) {
    deprecate_warn(
      "0.5.4",
      I(sprintf(
        "The %s argument(s) of the count loaders",
        paste0("`", names(supplied)[supplied], "`", collapse = ", ")
      )),
      details = paste(
        "The CSR to CSC conversion is bounded by `csc_mem_gb` instead.",
        "The old arguments are ignored."
      ),
      always = TRUE
    )
  }

  if (.verbose) {
    message(" Converting the cell-based data into the gene-based format.")
  }
  rust_con$generate_gene_based_data(max_mem_gb = csc_mem_gb, verbose = .verbose)

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
#' @param csc_mem_gb Optional numeric. Memory in GB for the buffers of the
#' cell-to-gene (CSR to CSC) conversion, at 10 bytes per non-zero. `NULL`
#' (default) converts in one pass and holds the whole matrix. Set a cap for
#' large data sets; every extra phase re-reads the cell file once.
#' @param streaming,batch_size,max_genes_in_memory,cell_batch_size Replaced by
#' `csc_mem_gb` and ignored. `r lifecycle::badge("deprecated")`
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @returns It will populate the files on disk and return the class with updated
#' shape information.
#'
#' @export
#'
#' @examplesIf requireNamespace("Seurat", quietly = TRUE)
#' \donttest{
#' # a Seurat object holds genes x cells, the loader transposes it
#' data <- generate_single_cell_test_data(
#'   syn_data_params = params_sc_synthetic_data(n_cells = 200L, n_genes = 40L)
#' )
#' seurat_obj <- Seurat::CreateSeuratObject(
#'   counts = Matrix::t(data$counts),
#'   meta.data = data.frame(data$obs, row.names = data$obs$cell_id)
#' )
#' dir_data <- tempfile("sc_seurat")
#' dir.create(dir_data, recursive = TRUE)
#' sc <- load_seurat(
#'   object = SingleCells(dir_data = dir_data),
#'   seurat = seurat_obj,
#'   sc_qc_param = params_sc_min_quality(
#'     min_unique_genes = 5L,
#'     min_lib_size = 25L,
#'     min_cells = 5L
#'   ),
#'   .verbose = FALSE
#' )
#' dim(sc)
#'
#' unlink(dir_data, recursive = TRUE, force = TRUE)
#' }
load_seurat <- S7::new_generic(
  name = "load_seurat",
  dispatch_args = "object",
  fun = function(
    object,
    seurat,
    sc_qc_param = params_sc_min_quality(),
    csc_mem_gb = NULL,
    streaming = deprecated(),
    batch_size = deprecated(),
    max_genes_in_memory = deprecated(),
    cell_batch_size = deprecated(),
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method load_seurat SingleCells
S7::method(load_seurat, SingleCells) <- function(
  object,
  seurat,
  sc_qc_param = params_sc_min_quality(),
  csc_mem_gb = NULL,
  streaming = deprecated(),
  batch_size = deprecated(),
  max_genes_in_memory = deprecated(),
  cell_batch_size = deprecated(),
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::assertClass(seurat, "Seurat")
  assertScMinQCParams(sc_qc_param)
  checkmate::qassert(csc_mem_gb, c("0", "N1(0,)"))
  checkmate::qassert(.verbose, "B1")

  if (.verbose) {
    message("Pulling the raw counts out of the Seurat object.")
  }

  counts <- get_seurat_counts(seurat)

  if (.verbose) {
    message("Pulling the obs and var data out of the object")
  }

  obs_dt <- data.table::as.data.table(
    seurat@meta.data,
    keep.rownames = "barcode"
  )

  # the first obs column becomes cell_id downstream, so an existing one would
  # collide and get a suffix from make.unique()
  if ("cell_id" %in% names(obs_dt)) {
    obs_dt[, cell_id := NULL]
  }

  var_dt <- data.table::data.table(gene_id = rownames(seurat))

  load_r_data(
    object = object,
    counts = counts,
    obs = obs_dt,
    var = var_dt,
    sc_qc_param = sc_qc_param,
    csc_mem_gb = csc_mem_gb,
    streaming = streaming,
    batch_size = batch_size,
    max_genes_in_memory = max_genes_in_memory,
    cell_batch_size = cell_batch_size,
    .verbose = .verbose
  )
}

## single cell experiment ------------------------------------------------------

#' Extract the obs and var tables from a SingleCellExperiment
#'
#' @description
#' `colData` becomes obs and `rowData` becomes var, with the identifiers put
#' first because the first column of each is what becomes `cell_id` and
#' `gene_id` downstream.
#'
#' Non-atomic columns are dropped. `colData` is a `DataFrame` and can carry
#' nested `DataFrame` or list columns, which have nowhere to go in the DuckDB.
#'
#' @param sce `SingleCellExperiment` class.
#'
#' @returns A list with `obs` and `var` as data.tables.
#'
#' @keywords internal
.sce_obs_var <- function(sce) {
  flatten <- function(df, ids, id_label, reserved, axis, n) {
    atomic <- vapply(
      seq_len(ncol(df)),
      \(i) is.atomic(df[[i]]) || is.factor(df[[i]]),
      logical(1)
    )
    if (any(!atomic)) {
      warning(sprintf(
        "Dropping %i non-atomic %s column(s): %s.",
        sum(!atomic),
        axis,
        paste(colnames(df)[!atomic], collapse = ", ")
      ))
    }

    kept <- if (any(atomic)) {
      data.table::as.data.table(as.list(df[, atomic, drop = FALSE]))
    } else {
      data.table::data.table()
    }

    # Bioconductor objects quite happily carry the identifiers in the metadata
    # and leave the dimnames empty, so fall back to the first column rather
    # than handing back a synthetic index nobody can join on
    if (is.null(ids)) {
      if (ncol(kept) > 0L) {
        ids <- as.character(kept[[1L]])
        warning(sprintf(
          "No %s names on the object. Using '%s' from the metadata instead.",
          axis,
          names(kept)[1L]
        ))
      } else {
        ids <- sprintf("%s_%i", axis, seq_len(n))
        warning(sprintf(
          "No %s names and no metadata to fall back on. Generating them.",
          axis
        ))
      }
    }

    # Identifiers have to be present and unique: they are the join key for
    # every lookup on this side. Published objects routinely fail both, e.g.
    # the ageing thymus data carries 336 genes with no annotation at all, so
    # repair rather than refuse, and say exactly what was repaired.
    ids <- as.character(ids)
    missing_ids <- is.na(ids) | !nzchar(ids)
    if (any(missing_ids)) {
      warning(sprintf(
        "%i %s(s) have no identifier. Generating one for each.",
        sum(missing_ids),
        axis
      ))
      ids[missing_ids] <- sprintf("%s_%i", axis, which(missing_ids))
    }
    if (anyDuplicated(ids)) {
      warning(sprintf(
        "%i %s identifier(s) are duplicated. Making them unique.",
        sum(duplicated(ids)),
        axis
      ))
      ids <- make.unique(ids)
    }

    # the leading column becomes `reserved` downstream, so anything that snake
    # cases to the same name collides and picks up a make.unique() suffix
    if (ncol(kept) > 0L) {
      clash <- to_snake_case(names(kept)) == reserved
      if (any(clash)) {
        kept[, (names(kept)[clash]) := NULL]
      }
    }

    out <- data.table::data.table(id = ids)
    data.table::setnames(out, "id", id_label)

    if (ncol(kept) > 0L) {
      out <- cbind(out, kept)
    }

    out
  }

  list(
    obs = flatten(
      SummarizedExperiment::colData(sce),
      colnames(sce),
      id_label = "barcode",
      reserved = "cell_id",
      axis = "cell",
      n = ncol(sce)
    ),
    var = flatten(
      SummarizedExperiment::rowData(sce),
      rownames(sce),
      id_label = "gene_id",
      reserved = "gene_id",
      axis = "gene",
      n = nrow(sce)
    )
  )
}

#' Load in data from a `SingleCellExperiment`
#'
#' @description
#' Brings a Bioconductor `SingleCellExperiment` into a `SingleCells` object.
#' `colData` becomes the obs table, `rowData` becomes the var table, and the
#' chosen assay goes through the same Rust quality control and normalisation
#' every other loader uses.
#'
#' The assay has to hold raw counts. Plenty of objects in the wild ship only
#' `logcounts`, and a negative binomial cannot model those, so pick the right
#' one rather than letting the default find whatever is there.
#'
#' `reducedDims` and `altExps` are not carried over. Run the embedding on this
#' side, and use [bixverse::SingleCellsMultiModal()] for ADT.
#'
#' @param object `SingleCells` class.
#' @param sce `SingleCellExperiment` class you want to transform.
#' @param assay_name String. Which assay holds the raw counts. Defaults to
#' `"counts"`.
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
#' @param csc_mem_gb Optional numeric. Memory in GB for the buffers of the
#' cell-to-gene (CSR to CSC) conversion, at 10 bytes per non-zero. `NULL`
#' (default) converts in one pass and holds the whole matrix. Set a cap for
#' large data sets; every extra phase re-reads the cell file once.
#' @param streaming,batch_size,max_genes_in_memory,cell_batch_size Replaced by
#' `csc_mem_gb` and ignored. `r lifecycle::badge("deprecated")`
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @returns It will populate the files on disk and return the class with updated
#' shape information.
#'
#' @export
#'
#' @examplesIf requireNamespace("SingleCellExperiment", quietly = TRUE)
#' \donttest{
#' # colData becomes obs, rowData becomes var, the counts assay gets normalised
#' data <- generate_single_cell_test_data(
#'   syn_data_params = params_sc_synthetic_data(n_cells = 200L, n_genes = 40L)
#' )
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'   assays = list(counts = as(Matrix::t(data$counts), "CsparseMatrix")),
#'   colData = data.frame(data$obs, row.names = data$obs$cell_id),
#'   rowData = data.frame(data$var, row.names = data$var$gene_id)
#' )
#' dir_data <- tempfile("sc_sce")
#' dir.create(dir_data, recursive = TRUE)
#' sc <- load_sce(
#'   object = SingleCells(dir_data = dir_data),
#'   sce = sce,
#'   sc_qc_param = params_sc_min_quality(
#'     min_unique_genes = 5L,
#'     min_lib_size = 25L,
#'     min_cells = 5L
#'   ),
#'   .verbose = FALSE
#' )
#' dim(sc)
#'
#' unlink(dir_data, recursive = TRUE, force = TRUE)
#' }
load_sce <- S7::new_generic(
  name = "load_sce",
  dispatch_args = "object",
  fun = function(
    object,
    sce,
    assay_name = "counts",
    sc_qc_param = params_sc_min_quality(),
    csc_mem_gb = NULL,
    streaming = deprecated(),
    batch_size = deprecated(),
    max_genes_in_memory = deprecated(),
    cell_batch_size = deprecated(),
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method load_sce SingleCells
S7::method(load_sce, SingleCells) <- function(
  object,
  sce,
  assay_name = "counts",
  sc_qc_param = params_sc_min_quality(),
  csc_mem_gb = NULL,
  streaming = deprecated(),
  batch_size = deprecated(),
  max_genes_in_memory = deprecated(),
  cell_batch_size = deprecated(),
  .verbose = TRUE
) {
  # checks
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop(
      paste(
        "Package 'SingleCellExperiment' required.",
        "Install with: BiocManager::install('SingleCellExperiment')"
      )
    )
  }

  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::assertClass(sce, "SingleCellExperiment")
  checkmate::qassert(assay_name, "S1")
  assertScMinQCParams(sc_qc_param)
  checkmate::qassert(csc_mem_gb, c("0", "N1(0,)"))
  checkmate::qassert(.verbose, "B1")

  available <- SummarizedExperiment::assayNames(sce)
  if (!assay_name %in% available) {
    stop(sprintf(
      "Assay '%s' not found. The object has: %s.",
      assay_name,
      paste(available, collapse = ", ")
    ))
  }

  if (.verbose) {
    message("Pulling the raw counts out of the SingleCellExperiment.")
  }

  counts <- SummarizedExperiment::assay(sce, assay_name)

  # a DelayedArray or HDF5Matrix has no slots to reinterpret, and realising it
  # silently would pull the whole thing into memory behind the user's back
  if (!inherits(counts, "dgCMatrix")) {
    stop(sprintf(
      paste(
        "Assay '%s' is a %s, not a dgCMatrix.",
        "Coerce it first, e.g. as(assay(sce, '%s'), 'CsparseMatrix'),",
        "and be aware that realising a DelayedArray loads it into memory."
      ),
      assay_name,
      class(counts)[1],
      assay_name
    ))
  }

  if (length(counts@x) > 0 && min(counts@x) < 0) {
    stop(sprintf(
      paste(
        "Assay '%s' holds negative values, so it is not raw counts.",
        "Pass the assay holding the counts via `assay_name`."
      ),
      assay_name
    ))
  }

  counts <- .counts_to_cell_major(counts)

  if (.verbose) {
    message("Pulling the obs and var data out of the object")
  }

  tables <- .sce_obs_var(sce)

  load_r_data(
    object = object,
    counts = counts,
    obs = tables$obs,
    var = tables$var,
    sc_qc_param = sc_qc_param,
    csc_mem_gb = csc_mem_gb,
    streaming = streaming,
    batch_size = batch_size,
    max_genes_in_memory = max_genes_in_memory,
    cell_batch_size = cell_batch_size,
    .verbose = .verbose
  )
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
#' @param csc_mem_gb Optional numeric. Memory in GB for the buffers of the
#' cell-to-gene (CSR to CSC) conversion, at 10 bytes per non-zero. `NULL`
#' (default) converts in one pass and holds the whole matrix. Set a cap for
#' large data sets; every extra phase re-reads the cell file once.
#' @param streaming,batch_size,max_genes_in_memory,cell_batch_size Replaced by
#' `csc_mem_gb` and ignored. `r lifecycle::badge("deprecated")`
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @returns It will populate the files on disk and return the class with updated
#' shape information.
#'
#' @export
#'
#' @examples
#' # straight from a dgRMatrix in memory onto disk
#' data <- generate_single_cell_test_data(
#'   syn_data_params = params_sc_synthetic_data(n_cells = 200L, n_genes = 40L)
#' )
#' dir_data <- tempfile("sc_r_data")
#' dir.create(dir_data, recursive = TRUE)
#' sc <- load_r_data(
#'   object = SingleCells(dir_data = dir_data),
#'   counts = data$counts,
#'   obs = data$obs,
#'   var = data$var,
#'   sc_qc_param = params_sc_min_quality(
#'     min_unique_genes = 5L,
#'     min_lib_size = 25L,
#'     min_cells = 5L
#'   ),
#'   .verbose = FALSE
#' )
#' dim(sc)
#'
#' unlink(dir_data, recursive = TRUE, force = TRUE)
load_r_data <- S7::new_generic(
  name = "load_r_data",
  dispatch_args = "object",
  fun = function(
    object,
    counts,
    obs,
    var,
    sc_qc_param = params_sc_min_quality(),
    csc_mem_gb = NULL,
    streaming = deprecated(),
    batch_size = deprecated(),
    max_genes_in_memory = deprecated(),
    cell_batch_size = deprecated(),
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method load_r_data SingleCells
S7::method(load_r_data, SingleCells) <- function(
  object,
  counts,
  obs,
  var,
  sc_qc_param = params_sc_min_quality(),
  csc_mem_gb = NULL,
  streaming = deprecated(),
  batch_size = deprecated(),
  max_genes_in_memory = deprecated(),
  cell_batch_size = deprecated(),
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::assertClass(counts, "dgRMatrix")
  no_cells <- nrow(counts)
  no_genes <- ncol(counts)
  checkmate::assertDataTable(obs, nrows = no_cells)
  checkmate::assertDataTable(var, nrows = no_genes)
  checkmate::qassert(csc_mem_gb, c("0", "N1(0,)"))
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
    csc_mem_gb = csc_mem_gb,
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

#### general -------------------------------------------------------------------

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
#' @param raw_count_slot Where raw counts live. `"auto"` detects per file via
#' [detect_raw_count_slot()]; otherwise one of `"X"`, `"raw.X"`,
#' `"layers.counts"`.
#' @param h5ad_streaming Boolean. Stream the h5ad counts into the cell-based
#' binary in batches instead of materialising the filtered matrix first.
#' Recommended for large files. Defaults to `TRUE`.
#' @param csc_mem_gb Optional numeric. Memory in GB for the buffers of the
#' cell-to-gene (CSR to CSC) conversion, at 10 bytes per non-zero. `NULL`
#' (default) converts in one pass and holds the whole matrix. Set a cap for
#' large data sets; every extra phase re-reads the cell file once.
#' @param streaming,batch_size,max_genes_in_memory,cell_batch_size Replaced by
#' `csc_mem_gb` and ignored. `r lifecycle::badge("deprecated")`
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @returns It will populate the files on disk and return the class with updated
#' shape information.
#'
#' @export
#'
#' @examples
#' # round trip through a sparse h5ad file
#' data <- generate_single_cell_test_data(
#'   syn_data_params = params_sc_synthetic_data(n_cells = 200L, n_genes = 40L)
#' )
#' f_path <- tempfile(fileext = ".h5ad")
#' write_h5ad_sc(f_path, data$counts, data$obs, data$var, .verbose = FALSE)
#' dir_data <- tempfile("sc_h5ad")
#' dir.create(dir_data, recursive = TRUE)
#' sc <- load_h5ad(
#'   object = SingleCells(dir_data = dir_data),
#'   h5_path = f_path,
#'   sc_qc_param = params_sc_min_quality(
#'     min_unique_genes = 5L,
#'     min_lib_size = 25L,
#'     min_cells = 5L
#'   ),
#'   .verbose = FALSE
#' )
#' dim(sc)
#'
#' unlink(c(f_path, dir_data), recursive = TRUE, force = TRUE)
load_h5ad <- S7::new_generic(
  name = "load_h5ad",
  dispatch_args = "object",
  fun = function(
    object,
    h5_path,
    sc_qc_param = params_sc_min_quality(),
    raw_count_slot = c("auto", "X", "raw.X", "layers.counts"),
    cell_id_col = NULL,
    h5ad_streaming = TRUE,
    csc_mem_gb = NULL,
    streaming = deprecated(),
    batch_size = deprecated(),
    max_genes_in_memory = deprecated(),
    cell_batch_size = deprecated(),
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method load_h5ad SingleCells
#'
#' @export
#'
#' @importFrom zeallot %<-%
#' @importFrom magrittr %>%
S7::method(load_h5ad, SingleCells) <- function(
  object,
  h5_path,
  sc_qc_param = params_sc_min_quality(),
  raw_count_slot = c("auto", "X", "raw.X", "layers.counts"),
  cell_id_col = NULL,
  h5ad_streaming = TRUE,
  csc_mem_gb = NULL,
  streaming = deprecated(),
  batch_size = deprecated(),
  max_genes_in_memory = deprecated(),
  cell_batch_size = deprecated(),
  .verbose = TRUE
) {
  raw_count_slot <- match.arg(raw_count_slot)

  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  assertScMinQCParams(sc_qc_param)
  checkmate::qassert(h5ad_streaming, "B1")
  checkmate::qassert(csc_mem_gb, c("0", "N1(0,)"))
  checkmate::assertChoice(
    raw_count_slot,
    c("auto", "X", "raw.X", "layers.counts")
  )
  checkmate::qassert(.verbose, "B1")

  raw_count_slot <- if (raw_count_slot == "auto") {
    resolved_slot <- detect_raw_count_slot(h5_path)
    if (is.na(resolved_slot)) {
      stop(paste(
        "No raw count slot could be found in the object!",
        "Please validate the h5ad file"
      ))
    }
    resolved_slot
  } else {
    raw_count_slot
  }

  h5_path <- path.expand(h5_path)

  h5_meta <- get_h5ad_dimensions(f_path = h5_path)

  rust_con <- get_sc_rust_ptr(object)

  h5ad_to_file <- if (h5ad_streaming) {
    rust_con$h5ad_to_file_streaming
  } else {
    rust_con$h5ad_to_file
  }

  file_res <- h5ad_to_file(
    cs_type = h5_meta$type,
    h5_path = h5_path,
    no_cells = h5_meta$dims["obs"],
    no_genes = h5_meta$dims["var"],
    qc_params = sc_qc_param,
    slot = raw_count_slot,
    verbose = .verbose
  )

  .dispatch_gene_based_data(
    rust_con = rust_con,
    csc_mem_gb = csc_mem_gb,
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
#' @param csc_mem_gb Optional numeric. Memory in GB for the buffers of the
#' cell-to-gene (CSR to CSC) conversion, at 10 bytes per non-zero. `NULL`
#' (default) converts in one pass and holds the whole matrix. Set a cap for
#' large data sets; every extra phase re-reads the cell file once.
#' @param streaming,batch_size,max_genes_in_memory,cell_batch_size Replaced by
#' `csc_mem_gb` and ignored. `r lifecycle::badge("deprecated")`
#' @param .verbose Boolean.
#'
#' @returns It will populate the files on disk and return the class with updated
#' shape information.
#'
#' @export
#'
#' @examples
#' # an h5ad holding log1p(x / lib_size * 1e4), raw counts reconstructed on read
#' data <- generate_single_cell_test_data(
#'   syn_data_params = params_sc_synthetic_data(n_cells = 200L, n_genes = 40L)
#' )
#' lib_size <- Matrix::rowSums(data$counts)
#' norm_counts <- data$counts
#' norm_counts@x <- log1p(
#'   norm_counts@x / rep(lib_size, diff(norm_counts@p)) * 1e4
#' )
#' obs <- data.table::copy(data$obs)[, total_counts := lib_size]
#'
#' f_path <- tempfile(fileext = ".h5ad")
#' write_h5ad_sc(f_path, norm_counts, obs, data$var, .verbose = FALSE)
#' dir_data <- tempfile("sc_h5ad_norm")
#' dir.create(dir_data, recursive = TRUE)
#' sc <- load_h5ad_norm(
#'   object = SingleCells(dir_data = dir_data),
#'   h5_path = f_path,
#'   obs_lib_size_col = "total_counts",
#'   target_size = 1e4,
#'   sc_qc_param = params_sc_min_quality(
#'     min_unique_genes = 5L,
#'     min_lib_size = 25L,
#'     min_cells = 5L
#'   ),
#'   .verbose = FALSE
#' )
#' dim(sc)
#'
#' unlink(c(f_path, dir_data), recursive = TRUE, force = TRUE)
load_h5ad_norm <- S7::new_generic(
  name = "load_h5ad_norm",
  dispatch_args = "object",
  fun = function(
    object,
    h5_path,
    obs_lib_size_col,
    target_size,
    sc_qc_param = params_sc_min_quality(),
    cell_id_col = NULL,
    csc_mem_gb = NULL,
    streaming = deprecated(),
    batch_size = deprecated(),
    max_genes_in_memory = deprecated(),
    cell_batch_size = deprecated(),
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method load_h5ad_norm SingleCells
#'
#' @export
#'
#' @importFrom zeallot %<-%
#' @importFrom magrittr %>%
S7::method(load_h5ad_norm, SingleCells) <- function(
  object,
  h5_path,
  obs_lib_size_col,
  target_size,
  sc_qc_param = params_sc_min_quality(),
  cell_id_col = NULL,
  csc_mem_gb = NULL,
  streaming = deprecated(),
  batch_size = deprecated(),
  max_genes_in_memory = deprecated(),
  cell_batch_size = deprecated(),
  .verbose = TRUE
) {
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  assertScMinQCParams(sc_qc_param)
  checkmate::qassert(obs_lib_size_col, "S1")
  checkmate::qassert(target_size, "N1")
  checkmate::qassert(csc_mem_gb, c("0", "N1(0,)"))
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
    csc_mem_gb = csc_mem_gb,
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
#' Convenience alias for `load_h5ad(h5ad_streaming = TRUE)`. Kept for
#' backwards compatibility. Prefer calling [bixverse::load_h5ad()] directly.
#'
#' @param object `SingleCells` class.
#' @param h5_path File path to the h5ad object.
#' @param sc_qc_param List. Output of [bixverse::params_sc_min_quality()].
#' @param raw_count_slot Where raw counts live. `"auto"` detects per file via
#' [detect_raw_count_slot()]; otherwise one of `"X"`, `"raw.X"`,
#' `"layers.counts"`.
#' @param csc_mem_gb Optional numeric. Memory in GB for the buffers of the
#' cell-to-gene (CSR to CSC) conversion, at 10 bytes per non-zero. `NULL`
#' (default) converts in one pass and holds the whole matrix. Set a cap for
#' large data sets; every extra phase re-reads the cell file once.
#' @param streaming,batch_size,max_genes_in_memory,cell_batch_size Replaced by
#' `csc_mem_gb` and ignored. `r lifecycle::badge("deprecated")`
#' @param .verbose Boolean.
#'
#' @returns The class with updated shape information.
#'
#' @export
#'
#' @examples
#' # same as load_h5ad(h5ad_streaming = TRUE)
#' data <- generate_single_cell_test_data(
#'   syn_data_params = params_sc_synthetic_data(n_cells = 200L, n_genes = 40L)
#' )
#' f_path <- tempfile(fileext = ".h5ad")
#' write_h5ad_sc(f_path, data$counts, data$obs, data$var, .verbose = FALSE)
#' dir_data <- tempfile("sc_stream")
#' dir.create(dir_data, recursive = TRUE)
#' sc <- stream_h5ad(
#'   object = SingleCells(dir_data = dir_data),
#'   h5_path = f_path,
#'   sc_qc_param = params_sc_min_quality(
#'     min_unique_genes = 5L,
#'     min_lib_size = 25L,
#'     min_cells = 5L
#'   ),
#'   .verbose = FALSE
#' )
#' dim(sc)
#'
#' unlink(c(f_path, dir_data), recursive = TRUE, force = TRUE)
stream_h5ad <- S7::new_generic(
  name = "stream_h5ad",
  dispatch_args = "object",
  fun = function(
    object,
    h5_path,
    sc_qc_param = params_sc_min_quality(),
    raw_count_slot = c("auto", "X", "raw.X", "layers.counts"),
    csc_mem_gb = NULL,
    streaming = deprecated(),
    batch_size = deprecated(),
    max_genes_in_memory = deprecated(),
    cell_batch_size = deprecated(),
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
  raw_count_slot = c("auto", "X", "raw.X", "layers.counts"),
  csc_mem_gb = NULL,
  streaming = deprecated(),
  batch_size = deprecated(),
  max_genes_in_memory = deprecated(),
  cell_batch_size = deprecated(),
  .verbose = TRUE
) {
  load_h5ad(
    object = object,
    h5_path = h5_path,
    sc_qc_param = sc_qc_param,
    raw_count_slot = raw_count_slot,
    h5ad_streaming = TRUE,
    csc_mem_gb = csc_mem_gb,
    streaming = streaming,
    batch_size = batch_size,
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
#' @param csc_mem_gb Optional numeric. Memory in GB for the buffers of the
#' cell-to-gene (CSR to CSC) conversion, at 10 bytes per non-zero. `NULL`
#' (default) converts in one pass and holds the whole matrix. Set a cap for
#' large data sets; every extra phase re-reads the cell file once.
#' @param streaming,batch_size,max_genes_in_memory,cell_batch_size Replaced by
#' `csc_mem_gb` and ignored. `r lifecycle::badge("deprecated")`
#' @param .verbose Boolean.
#'
#' @returns The class with updated shape and populated DuckDB.
#'
#' @export
#'
#' @examples
#' # two files into one experiment, cells tagged by exp_id
#' files <- c(a = tempfile(fileext = ".h5ad"), b = tempfile(fileext = ".h5ad"))
#' for (i in seq_along(files)) {
#'   data <- generate_single_cell_test_data(
#'     syn_data_params = params_sc_synthetic_data(
#'       n_cells = 200L,
#'       n_genes = 40L
#'     ),
#'     seed = i
#'   )
#'   write_h5ad_sc(files[i], data$counts, data$obs, data$var, .verbose = FALSE)
#' }
#' tasks <- prescan_h5ad_files(h5_paths = files, .verbose = FALSE)
#' dir_data <- tempfile("sc_multi_h5ad")
#' dir.create(dir_data, recursive = TRUE)
#' sc <- load_multi_h5ad(
#'   object = SingleCells(dir_data = dir_data),
#'   prescan_result = tasks,
#'   sc_qc_param = params_sc_min_quality(
#'     min_unique_genes = 5L,
#'     min_lib_size = 25L,
#'     min_cells = 5L
#'   ),
#'   .verbose = FALSE
#' )
#' table(sc[["exp_id"]])
#'
#' unlink(c(files, dir_data), recursive = TRUE, force = TRUE)
load_multi_h5ad <- S7::new_generic(
  name = "load_multi_h5ad",
  dispatch_args = "object",
  fun = function(
    object,
    prescan_result,
    sc_qc_param = params_sc_min_quality(),
    cell_id_col = NULL,
    csc_mem_gb = NULL,
    streaming = deprecated(),
    batch_size = deprecated(),
    max_genes_in_memory = deprecated(),
    cell_batch_size = deprecated(),
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
  csc_mem_gb = NULL,
  streaming = deprecated(),
  batch_size = deprecated(),
  max_genes_in_memory = deprecated(),
  cell_batch_size = deprecated(),
  .verbose = TRUE
) {
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  assertScMinQCParams(sc_qc_param)
  checkmate::assertList(prescan_result)
  checkmate::assertTRUE(all(
    c("universe", "universe_size", "file_tasks") %in%
      names(prescan_result)
  ))
  checkmate::qassert(csc_mem_gb, c("0", "N1(0,)"))
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
    csc_mem_gb = csc_mem_gb,
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
  duckdb_con$populate_var_minimal(final_gene_names = final_gene_names)

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

#### export --------------------------------------------------------------------

#' Save a `SingleCells` object to h5ad
#'
#' @description
#' Writes the counts, the DuckDB obs/var tables and whatever is cached in
#' memory (PCA, embeddings, sNN graph) to a spec-compliant h5ad file, so that
#' the experiment can be handed to ScanPy or shared. The counts are streamed
#' cell batch by cell batch and never fully materialised in R.
#'
#' Only the cells that are currently kept are written, i.e. the export matches
#' what `object[]` and `object[[]]` return. Counts are stored as `float32`,
#' which is the ScanPy convention.
#'
#' A `SingleCellsMultiModal` object inherits this method and exports its RNA
#' modality; the ADT layer is not written.
#'
#' @param object `SingleCells` class.
#' @param h5_path File path to write the h5ad file to.
#' @param assay String. One of `c("raw", "norm")`. Which count assay to place
#' in `X`.
#' @param chunk_size Integer. Number of cells per streaming batch. Defaults to
#' `10000L`.
#' @param overwrite Boolean. Shall an existing file be overwritten. Defaults to
#' `TRUE`.
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @return Returns the path to the written file, invisibly.
#'
#' @export
save_h5ad <- S7::new_generic(
  name = "save_h5ad",
  dispatch_args = "object",
  fun = function(
    object,
    h5_path,
    assay = c("raw", "norm"),
    chunk_size = 10000L,
    overwrite = TRUE,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method save_h5ad SingleCells
#'
#' @export
S7::method(save_h5ad, SingleCells) <- function(
  object,
  h5_path,
  assay = c("raw", "norm"),
  chunk_size = 10000L,
  overwrite = TRUE,
  .verbose = TRUE
) {
  assay <- match.arg(assay)

  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::qassert(h5_path, "S1")
  checkmate::assertChoice(assay, c("raw", "norm"))
  checkmate::qassert(chunk_size, "I1")
  checkmate::assertTRUE(chunk_size > 0L)
  checkmate::qassert(overwrite, "B1")
  checkmate::qassert(.verbose, "B1")

  h5_path <- path.expand(h5_path)
  checkmate::assertPathForOutput(h5_path, overwrite = TRUE)

  if (file.exists(h5_path)) {
    if (!overwrite) {
      stop("The h5ad file already exists and overwrite = FALSE.")
    }
    file.remove(h5_path)
  }

  obs <- .h5ad_export_obs(object)
  var <- .h5ad_export_var(object)

  # `set_cells_to_keep()` keeps whatever order it was handed, while the obs
  # table comes back in cell index order; sorting here is what keeps the two
  # aligned. 0-based, as Rust wants it.
  cell_indices <- sort(get_cells_to_keep(object))

  # the DuckDB and the ScMap both track which cells are kept; a mismatch here
  # would silently misalign the obs table against the counts
  checkmate::assertTRUE(nrow(obs) == length(cell_indices))
  checkmate::assertTRUE(nrow(var) == S7::prop(object, "dims")[2L])

  cache <- .h5ad_export_cache(object = object, .verbose = .verbose)

  if (.verbose) {
    message(sprintf(
      "Writing %i cells x %i genes to h5ad.",
      nrow(obs),
      nrow(var)
    ))
  }

  rs_save_h5ad(
    f_path_cells = get_rust_count_cell_f_path(object),
    h5_path = h5_path,
    cell_indices = as.integer(cell_indices),
    norm = assay == "norm",
    obs_index = as.character(obs[["cell_id"]]),
    obs = .h5ad_columns(obs[, !"cell_id"]),
    var_index = as.character(var[["gene_id"]]),
    var = .h5ad_columns(var[, !"gene_id"]),
    obsm = cache$obsm,
    varm = cache$varm,
    obsp = cache$obsp,
    uns_json = cache$uns_json,
    chunk_size = chunk_size
  )

  invisible(h5_path)
}

##### export helpers -----------------------------------------------------------

#' Assemble the obs table for an h5ad export
#'
#' @param object `SingleCells` class.
#'
#' @return A data.table with `cell_id` first and the bixverse bookkeeping
#' columns dropped.
#'
#' @keywords internal
.h5ad_export_obs <- function(object) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))

  obs <- data.table::copy(get_sc_obs(object, filtered = TRUE))

  # the counts are streamed in cell index order, so the obs table has to be
  if ("cell_idx" %in% names(obs)) {
    data.table::setorderv(obs, "cell_idx")
  }

  drop_cols <- intersect(c("cell_idx", "to_keep"), names(obs))
  if (length(drop_cols) > 0L) {
    obs[, (drop_cols) := NULL]
  }

  data.table::setcolorder(obs, "cell_id")

  return(obs)
}

#' Assemble the var table for an h5ad export
#'
#' @param object `SingleCells` class.
#'
#' @return A data.table with `gene_id` first and the bixverse bookkeeping
#' columns dropped.
#'
#' @keywords internal
.h5ad_export_var <- function(object) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))

  var <- data.table::copy(get_sc_var(object))

  if ("gene_idx" %in% names(var)) {
    data.table::setorderv(var, "gene_idx")
    var[, "gene_idx" := NULL]
  }

  data.table::setcolorder(var, "gene_id")

  return(var)
}

#' Bring table columns into the types the h5ad writer encodes
#'
#' @description The Rust side maps factors to categoricals, characters to
#' string arrays, and doubles, integers and logicals to plain arrays. This
#' picks which of those each column becomes:
#'
#' - a character column with repeated values becomes a factor, which is what
#'   pandas would hold it as
#' - a column with nothing but missing values has no category to write and
#'   becomes the literal string `"NA"`
#' - integers and logicals with missing values are widened to double, the
#'   only one of the three with a representation for them (`NA_integer_` is
#'   `INT_MIN` on disk, a silently wrong number)
#' - anything else (dates, lists) is stringified rather than dropped
#'
#' @param dt data.table. The columns to write, without the index.
#'
#' @return A named list of factor, character, double, integer or logical
#' vectors.
#'
#' @keywords internal
.h5ad_columns <- function(dt) {
  # checks
  checkmate::assertDataTable(dt)

  purrr::map(as.list(dt), function(x) {
    if (is.character(x) && data.table::uniqueN(x) < length(x)) {
      x <- factor(x)
    }
    if (is.factor(x)) {
      if (nlevels(x) == 0L) as.character(x) else x
    } else if (is.character(x)) {
      x
    } else if ((is.integer(x) || is.logical(x)) && anyNA(x)) {
      as.numeric(x)
    } else if (is.numeric(x) || is.logical(x)) {
      x
    } else {
      as.character(x)
    }
  })
}

#' Collect the cached embeddings and graphs of a `SingleCells` object
#'
#' @description PCA factors go to `obsm/X_pca` and the loadings to
#' `varm/PCs`, every other embedding to `obsm/X_<name>`, and the sNN graph to
#' `obsp/connectivities` together with the `uns/neighbors` block ScanPy looks
#' for. Anything that is not cached, or that no longer matches the cells being
#' written, is skipped.
#'
#' @param object `SingleCells` class.
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @return A list with `obsm` and `varm` (named lists of double matrices),
#' `obsp` (named list of CSR matrices as `indptr`, `indices`, `data`) and
#' `uns_json` (string).
#'
#' @keywords internal
.h5ad_export_cache <- function(object, .verbose = TRUE) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::qassert(.verbose, "B1")

  sc_cache <- get_sc_cache(object)
  no_cells <- length(get_cells_to_keep(object))
  no_genes <- S7::prop(object, "dims")[2L]

  obsm <- list()
  varm <- list()
  obsp <- list()
  uns <- list()

  embeddings <- setdiff(get_available_embeddings(sc_cache), "")

  for (embd_name in embeddings) {
    embd <- .drop_stamp(get_embedding(sc_cache, embd_name = embd_name))
    if (!is.matrix(embd)) {
      next
    }
    # a cached embedding can predate the current cell filter, in which case
    # anndata would reject the file; skip it rather than write something that
    # cannot be read back
    if (nrow(embd) != no_cells) {
      warning(sprintf(
        paste(
          "Embedding '%s' has %i rows but %i cells are being written.",
          "Skipping it in the h5ad export."
        ),
        embd_name,
        nrow(embd),
        no_cells
      ))
      next
    }
    h5_name <- sprintf("X_%s", embd_name)
    if (.verbose) {
      message(sprintf(
        "Adding the %s embedding as obsm/%s.",
        embd_name,
        h5_name
      ))
    }
    storage.mode(embd) <- "double"
    obsm[[h5_name]] <- unname(embd)
  }

  loadings <- .drop_stamp(get_pca_loadings(sc_cache))
  if (is.matrix(loadings)) {
    # the loadings only cover the HVGs, so they are padded back out to every
    # gene; ScanPy expects varm to be aligned with var. `get_hvg()` is 0-based
    # for Rust.
    hvg <- get_hvg(object)
    if (length(hvg) == nrow(loadings)) {
      padded <- matrix(0, nrow = no_genes, ncol = ncol(loadings))
      padded[hvg + 1L, ] <- loadings
      varm[["PCs"]] <- padded
    } else {
      warning(sprintf(
        paste(
          "The cached PCA loadings cover %i genes but %i HVGs are set.",
          "Skipping varm/PCs in the h5ad export."
        ),
        nrow(loadings),
        length(hvg)
      ))
    }
  }

  snn_graph <- .drop_stamp(get_snn_graph(sc_cache))
  if (inherits(snn_graph, "igraph") && igraph::vcount(snn_graph) == no_cells) {
    if (.verbose) {
      message("Adding the sNN graph as obsp/connectivities.")
    }
    attr_name <- if ("weight" %in% igraph::edge_attr_names(snn_graph)) {
      "weight"
    } else {
      NULL
    }
    adj <- igraph::as_adjacency_matrix(
      snn_graph,
      attr = attr_name,
      sparse = TRUE
    )
    adj <- methods::as(
      methods::as(methods::as(adj, "dMatrix"), "generalMatrix"),
      "RsparseMatrix"
    )
    obsp[["connectivities"]] <- list(
      indptr = adj@p,
      indices = adj@j,
      data = adj@x
    )
    uns[["neighbors"]] <- list(connectivities_key = "connectivities")
  }

  list(
    obsm = obsm,
    varm = varm,
    obsp = obsp,
    uns_json = as.character(jsonlite::toJSON(uns, auto_unbox = TRUE))
  )
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
#' @param csc_mem_gb Optional numeric. Memory in GB for the buffers of the
#' cell-to-gene (CSR to CSC) conversion, at 10 bytes per non-zero. `NULL`
#' (default) converts in one pass and holds the whole matrix. Set a cap for
#' large data sets; every extra phase re-reads the cell file once.
#' @param streaming,batch_size,max_genes_in_memory,cell_batch_size Replaced by
#' `csc_mem_gb` and ignored. `r lifecycle::badge("deprecated")`
#' @param .verbose Boolean.
#'
#' @returns The class with updated shape information.
#'
#' @export
#'
#' @examples
#' # read back a CellRanger style .mtx trio
#' data <- generate_single_cell_test_data(
#'   syn_data_params = params_sc_synthetic_data(n_cells = 200L, n_genes = 40L)
#' )
#' dir_src <- tempfile("cellranger")
#' dir.create(dir_src, recursive = TRUE)
#' write_cellranger_output(
#'   dir_src, data$counts, data$obs, data$var,
#'   rows = "cells", format_type = "csv", .verbose = FALSE
#' )
#'
#' dir_data <- tempfile("sc_mtx")
#' dir.create(dir_data, recursive = TRUE)
#' sc <- load_mtx(
#'   object = SingleCells(dir_data = dir_data),
#'   sc_mtx_io_param = params_sc_mtx_io(
#'     path_mtx = file.path(dir_src, "matrix.mtx"),
#'     path_obs = file.path(dir_src, "barcodes.csv"),
#'     path_var = file.path(dir_src, "features.csv"),
#'     cells_as_rows = TRUE,
#'     has_hdr = TRUE
#'   ),
#'   sc_qc_param = params_sc_min_quality(
#'     min_unique_genes = 5L,
#'     min_lib_size = 25L,
#'     min_cells = 5L
#'   ),
#'   .verbose = FALSE
#' )
#' dim(sc)
#'
#' unlink(c(dir_src, dir_data), recursive = TRUE, force = TRUE)
load_mtx <- S7::new_generic(
  name = "load_mtx",
  dispatch_args = "object",
  fun = function(
    object,
    sc_mtx_io_param = params_sc_mtx_io(),
    sc_qc_param = params_sc_min_quality(),
    mtx_streaming = TRUE,
    csc_mem_gb = NULL,
    streaming = deprecated(),
    batch_size = deprecated(),
    max_genes_in_memory = deprecated(),
    cell_batch_size = deprecated(),
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method load_mtx SingleCells
#'
#' @export
#'
#' @importFrom zeallot %<-%
#' @importFrom magrittr %>%
S7::method(load_mtx, SingleCells) <- function(
  object,
  sc_mtx_io_param = params_sc_mtx_io(),
  sc_qc_param = params_sc_min_quality(),
  mtx_streaming = TRUE,
  csc_mem_gb = NULL,
  streaming = deprecated(),
  batch_size = deprecated(),
  max_genes_in_memory = deprecated(),
  cell_batch_size = deprecated(),
  .verbose = TRUE
) {
  checkmate::assertClass(object, "bixverse::SingleCells")
  assertScMtxIOParams(sc_mtx_io_param)
  assertScMinQCParams(sc_qc_param)
  checkmate::qassert(mtx_streaming, "B1")
  checkmate::qassert(csc_mem_gb, c("0", "N1(0,)"))
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
    csc_mem_gb = csc_mem_gb,
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
#' @param csc_mem_gb Optional numeric. Memory in GB for the buffers of the
#' cell-to-gene (CSR to CSC) conversion, at 10 bytes per non-zero. `NULL`
#' (default) converts in one pass and holds the whole matrix. Set a cap for
#' large data sets; every extra phase re-reads the cell file once.
#' @param streaming,batch_size,max_genes_in_memory,cell_batch_size Replaced by
#' `csc_mem_gb` and ignored. `r lifecycle::badge("deprecated")`
#' @param .verbose Boolean.
#'
#' @returns The class with updated shape and populated DuckDB.
#'
#' @export
#'
#' @examples
#' # two CellRanger directories into one experiment
#' dirs <- c(tempfile("cr_a"), tempfile("cr_b"))
#' for (i in seq_along(dirs)) {
#'   dir.create(dirs[i], recursive = TRUE)
#'   data <- generate_single_cell_test_data(
#'     syn_data_params = params_sc_synthetic_data(
#'       n_cells = 200L,
#'       n_genes = 40L
#'     ),
#'     seed = i
#'   )
#'   write_cellranger_output(
#'     dirs[i], data$counts, data$obs, data$var,
#'     rows = "cells", format_type = "csv", .verbose = FALSE
#'   )
#' }
#' scan_res <- prescan_mtx_dirs(
#'   dirs = dirs,
#'   exp_ids = c("a", "b"),
#'   cells_as_rows = TRUE,
#'   has_hdr = TRUE,
#'   .verbose = FALSE
#' )
#'
#' dir_data <- tempfile("sc_multi_mtx")
#' dir.create(dir_data, recursive = TRUE)
#' sc <- load_multi_mtx(
#'   object = SingleCells(dir_data = dir_data),
#'   prescan_result = scan_res,
#'   sc_qc_param = params_sc_min_quality(
#'     min_unique_genes = 5L,
#'     min_lib_size = 25L,
#'     min_cells = 5L
#'   ),
#'   .verbose = FALSE
#' )
#' dim(sc)
#'
#' unlink(
#'   c(dirs, dir_data, scan_res$temp_files),
#'   recursive = TRUE,
#'   force = TRUE
#' )
load_multi_mtx <- S7::new_generic(
  name = "load_multi_mtx",
  dispatch_args = "object",
  fun = function(
    object,
    prescan_result,
    sc_qc_param = params_sc_min_quality(),
    csc_mem_gb = NULL,
    streaming = deprecated(),
    batch_size = deprecated(),
    max_genes_in_memory = deprecated(),
    cell_batch_size = deprecated(),
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
  csc_mem_gb = NULL,
  streaming = deprecated(),
  batch_size = deprecated(),
  max_genes_in_memory = deprecated(),
  cell_batch_size = deprecated(),
  .verbose = TRUE
) {
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  assertScMinQCParams(sc_qc_param)
  checkmate::assertList(prescan_result)
  checkmate::assertTRUE(all(
    c("universe", "universe_size", "file_tasks", "temp_files") %in%
      names(prescan_result)
  ))
  checkmate::qassert(csc_mem_gb, c("0", "N1(0,)"))
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
    csc_mem_gb = csc_mem_gb,
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
  duckdb_con$populate_var_minimal(final_gene_names = final_gene_names)

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

### h5 10x outputs -------------------------------------------------------------

#' Load in a 10x CellRanger h5 file to `SingleCells`
#'
#' @description
#' Loads the gene-expression modality from a CellRanger v2/v3 h5 file. The
#' counts go into the Rust-binarised format (with log normalisation applied on
#' read) and the barcodes/features into the DuckDB. Non-gene modalities (e.g.
#' Antibody Capture) are filtered out via `feature_type`.
#'
#' @param object `SingleCells` class.
#' @param h5_path File path to the 10x h5 file.
#' @param sc_qc_param List. Output of [bixverse::params_sc_min_quality()].
#' @param feature_type String. Modality to keep. Defaults to
#' `"Gene Expression"`. Ignored for v2 (single modality).
#' @param csc_mem_gb Optional numeric. Memory in GB for the buffers of the
#' cell-to-gene (CSR to CSC) conversion, at 10 bytes per non-zero. `NULL`
#' (default) converts in one pass and holds the whole matrix. Set a cap for
#' large data sets; every extra phase re-reads the cell file once.
#' @param streaming,batch_size,max_genes_in_memory,cell_batch_size Replaced by
#' `csc_mem_gb` and ignored. `r lifecycle::badge("deprecated")`
#' @param .verbose Boolean.
#'
#' @returns The class with updated shape information.
#'
#' @export
#'
#' @examples
#' # read back a CellRanger v3 h5, gene expression only
#' data <- generate_single_cell_test_data(
#'   syn_data_params = params_sc_synthetic_data(n_cells = 200L, n_genes = 40L)
#' )
#' f_path <- tempfile(fileext = ".h5")
#' write_tenx_h5_sc(
#'   f_path = f_path,
#'   counts = data$counts,
#'   barcodes = data$obs$cell_id,
#'   features = data.table::data.table(
#'     id = data$var$gene_id,
#'     name = data$var$ensembl_id,
#'     feature_type = "Gene Expression"
#'   )
#' )
#'
#' dir_data <- tempfile("sc_tenx")
#' dir.create(dir_data, recursive = TRUE)
#' sc <- load_tenx_h5(
#'   object = SingleCells(dir_data = dir_data),
#'   h5_path = f_path,
#'   sc_qc_param = params_sc_min_quality(
#'     min_unique_genes = 5L,
#'     min_lib_size = 25L,
#'     min_cells = 5L
#'   ),
#'   .verbose = FALSE
#' )
#' dim(sc)
#'
#' unlink(c(f_path, dir_data), recursive = TRUE, force = TRUE)
load_tenx_h5 <- S7::new_generic(
  name = "load_tenx_h5",
  dispatch_args = "object",
  fun = function(
    object,
    h5_path,
    sc_qc_param = params_sc_min_quality(),
    feature_type = "Gene Expression",
    csc_mem_gb = NULL,
    streaming = deprecated(),
    batch_size = deprecated(),
    max_genes_in_memory = deprecated(),
    cell_batch_size = deprecated(),
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method load_tenx_h5 SingleCells
#'
#' @export
S7::method(load_tenx_h5, SingleCells) <- function(
  object,
  h5_path,
  sc_qc_param = params_sc_min_quality(),
  feature_type = "Gene Expression",
  csc_mem_gb = NULL,
  streaming = deprecated(),
  batch_size = deprecated(),
  max_genes_in_memory = deprecated(),
  cell_batch_size = deprecated(),
  .verbose = TRUE
) {
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  assertScMinQCParams(sc_qc_param)
  checkmate::qassert(feature_type, c("S1", "0"))
  checkmate::qassert(csc_mem_gb, c("0", "N1(0,)"))
  checkmate::qassert(.verbose, "B1")

  h5_path <- path.expand(h5_path)
  meta <- get_tenx_h5_metadata(h5_path)

  rust_con <- get_sc_rust_ptr(object)

  file_res <- rust_con$tenx_h5_to_file_streaming(
    h5_path = h5_path,
    version = meta$version,
    no_cells = meta$n_cells,
    no_genes = meta$n_genes,
    qc_params = sc_qc_param,
    feature_type = feature_type,
    verbose = .verbose
  )

  .dispatch_gene_based_data(
    rust_con = rust_con,
    csc_mem_gb = csc_mem_gb,
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
    message("Loading barcodes from 10x h5 into the DuckDB.")
  }
  duckdb_con$populate_obs_from_tenx_h5(
    h5_path = h5_path,
    version = meta$version,
    filter = as.integer(file_res$cell_indices + 1)
  )

  if (.verbose) {
    message("Loading features from 10x h5 into the DuckDB.")
  }
  duckdb_con$populate_vars_from_tenx_h5(
    h5_path = h5_path,
    version = meta$version,
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

### multiple h5 10x outputs ----------------------------------------------------

#' Load multiple 10x CellRanger h5 files into a single `SingleCells`
#'
#' @description
#' Takes the result of [bixverse::prescan_tenx_h5_files()] and loads all
#' inputs into a single experiment with global gene QC and sequential cell
#' indexing. The feature space is determined by the prescan
#' (intersection or union of gene ids).
#'
#' @param object `SingleCells` class.
#' @param prescan_result Output of [bixverse::prescan_tenx_h5_files()].
#' @param sc_qc_param List. Output of [bixverse::params_sc_min_quality()].
#' @param csc_mem_gb Optional numeric. Memory in GB for the buffers of the
#' cell-to-gene (CSR to CSC) conversion, at 10 bytes per non-zero. `NULL`
#' (default) converts in one pass and holds the whole matrix. Set a cap for
#' large data sets; every extra phase re-reads the cell file once.
#' @param streaming,batch_size,max_genes_in_memory,cell_batch_size Replaced by
#' `csc_mem_gb` and ignored. `r lifecycle::badge("deprecated")`
#' @param .verbose Boolean.
#'
#' @returns The class with updated shape and populated DuckDB.
#'
#' @export
#'
#' @examples
#' # two 10x h5 files into one experiment
#' files <- c(a = tempfile(fileext = ".h5"), b = tempfile(fileext = ".h5"))
#' for (i in seq_along(files)) {
#'   data <- generate_single_cell_test_data(
#'     syn_data_params = params_sc_synthetic_data(
#'       n_cells = 200L,
#'       n_genes = 40L
#'     ),
#'     seed = i
#'   )
#'   write_tenx_h5_sc(
#'     f_path = files[i],
#'     counts = data$counts,
#'     barcodes = data$obs$cell_id,
#'     features = data.table::data.table(
#'       id = data$var$gene_id,
#'       name = data$var$ensembl_id,
#'       feature_type = "Gene Expression"
#'     )
#'   )
#' }
#' scan_res <- prescan_tenx_h5_files(h5_paths = files, .verbose = FALSE)
#'
#' dir_data <- tempfile("sc_multi_tenx")
#' dir.create(dir_data, recursive = TRUE)
#' sc <- load_multi_tenx_h5(
#'   object = SingleCells(dir_data = dir_data),
#'   prescan_result = scan_res,
#'   sc_qc_param = params_sc_min_quality(
#'     min_unique_genes = 5L,
#'     min_lib_size = 25L,
#'     min_cells = 5L
#'   ),
#'   .verbose = FALSE
#' )
#' dim(sc)
#'
#' unlink(c(files, dir_data), recursive = TRUE, force = TRUE)
load_multi_tenx_h5 <- S7::new_generic(
  name = "load_multi_tenx_h5",
  dispatch_args = "object",
  fun = function(
    object,
    prescan_result,
    sc_qc_param = params_sc_min_quality(),
    csc_mem_gb = NULL,
    streaming = deprecated(),
    batch_size = deprecated(),
    max_genes_in_memory = deprecated(),
    cell_batch_size = deprecated(),
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method load_multi_tenx_h5 SingleCells
#'
#' @export
S7::method(load_multi_tenx_h5, SingleCells) <- function(
  object,
  prescan_result,
  sc_qc_param = params_sc_min_quality(),
  csc_mem_gb = NULL,
  streaming = deprecated(),
  batch_size = deprecated(),
  max_genes_in_memory = deprecated(),
  cell_batch_size = deprecated(),
  .verbose = TRUE
) {
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  assertScMinQCParams(sc_qc_param)
  checkmate::assertList(prescan_result)
  checkmate::assertTRUE(all(
    c("universe", "universe_size", "file_tasks") %in% names(prescan_result)
  ))
  checkmate::qassert(csc_mem_gb, c("0", "N1(0,)"))
  checkmate::qassert(.verbose, "B1")

  rust_con <- get_sc_rust_ptr(object)

  file_res <- rust_con$multi_tenx_h5_to_file(
    file_tasks = prescan_result$file_tasks,
    universe_size = as.integer(prescan_result$universe_size),
    qc_params = sc_qc_param,
    verbose = .verbose
  )

  .dispatch_gene_based_data(
    rust_con = rust_con,
    csc_mem_gb = csc_mem_gb,
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
    message("Loading barcodes from 10x h5 files into DuckDB.")
  }

  per_file_info <- lapply(file_res$per_file, function(f) {
    task <- prescan_result$file_tasks[[f$exp_id]]
    list(
      h5_path = task$h5_path,
      version = task$version,
      exp_id = f$exp_id,
      cell_filter = as.integer(f$cell_indices + 1L)
    )
  })

  duckdb_con$populate_obs_from_multi_tenx_h5(per_file_info = per_file_info)

  if (.verbose) {
    message("Loading features into DuckDB.")
  }

  final_gene_names <- prescan_result$universe[file_res$global_gene_indices + 1L]
  duckdb_con$populate_var_minimal(final_gene_names = final_gene_names)

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

# generic in base_generics_sc.R

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

# generic in base_generics_sc.R

#' @method load_existing SingleCells
#'
#' @export
#'
#' @importFrom zeallot %<-%
#' @importFrom magrittr %>%
S7::method(load_existing, SingleCells) <- function(object, .verbose = TRUE) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCells))
  checkmate::qassert(.verbose, "B1")

  dir_data <- S7::prop(object, "dir_data")

  checkmate::assertFileExists(file.path(dir_data, "counts_cells.bin"))
  checkmate::assertFileExists(file.path(dir_data, "counts_genes.bin"))
  checkmate::assertFileExists(file.path(dir_data, "sc_duckdb.db"))

  # function body
  rust_con <- get_sc_rust_ptr(object)
  rust_con$set_from_file()

  duckdb_con <- get_sc_duckdb(object)

  if (any(c("memory.qs2", "memory.rds") %in% list.files(dir_data))) {
    if (.verbose) {
      message(paste(
        "Found stored data from save_sc_exp_to_disk().",
        "Loading that one into the object."
      ))
    }

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
    # objects saved before a slot was added come back without it
    S7::prop(object, "sc_map") <- .migrate_sc_map(saved_data$sc_map)
    S7::prop(object, "sc_cache") <- .migrate_sc_cache(saved_data$sc_cache)

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

# multi modal i/o --------------------------------------------------------------

## methods ---------------------------------------------------------------------

### save to disk ---------------------------------------------------------------

#' @method save_sc_exp_to_disk SingleCellsMultiModal
#'
#' @export
S7::method(save_sc_exp_to_disk, SingleCellsMultiModal) <- function(
  object,
  type = c("qs2", "rds")
) {
  type <- match.arg(type)
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCellsMultiModal))
  checkmate::assertChoice(type, c("qs2", "rds"))

  # pull the data out from the class
  sc_map <- get_sc_map(object)
  sc_cache <- get_sc_cache(object)
  adt_cache <- S7::prop(object, "adt_cache")
  atac_cache <- S7::prop(object, "atac_cache")
  adt_counts <- S7::prop(object, "adt_counts")
  other_data <- S7::prop(object, "other_data")
  dir <- S7::prop(object, "dir_data")

  to_save <- list(
    sc_map = sc_map,
    sc_cache = sc_cache,
    adt_cache = adt_cache,
    atac_cache = atac_cache,
    adt_counts = adt_counts,
    other_data = other_data
  )

  if (type == "qs2") {
    if (!requireNamespace("qs2", quietly = TRUE)) {
      stop("Package 'qs2' is required to use qs2 format. Please install it.")
    }
    qs2::qs_save(to_save, file = file.path(dir, "memory.qs2"))
  } else if (type == "rds") {
    saveRDS(to_save, file = file.path(dir, "memory.rds"))
  }
}

### load from disk -------------------------------------------------------------

#' @method load_existing SingleCellsMultiModal
#'
#' @export
#'
#' @importFrom zeallot %<-%
#' @importFrom magrittr %>%
S7::method(load_existing, SingleCellsMultiModal) <- function(
  object,
  .verbose = TRUE
) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, SingleCellsMultiModal))
  checkmate::qassert(.verbose, "B1")

  dir_data <- S7::prop(object, "dir_data")

  checkmate::assertFileExists(file.path(dir_data, "counts_cells.bin"))
  checkmate::assertFileExists(file.path(dir_data, "counts_genes.bin"))
  checkmate::assertFileExists(file.path(dir_data, "sc_duckdb.db"))

  # function body
  rust_con <- get_sc_rust_ptr(object)
  rust_con$set_from_file()

  duckdb_con <- get_sc_duckdb(object)

  if (any(c("memory.qs2", "memory.rds") %in% list.files(dir_data))) {
    if (.verbose) {
      message(paste(
        "Found stored data from save_sc_exp_to_disk().",
        "Loading that one into the object."
      ))
    }

    saved_data <- if ("memory.qs2" %in% list.files(dir_data)) {
      if (!requireNamespace("qs2", quietly = TRUE)) {
        stop("Package 'qs2' is required to use qs2 format. Please install it.")
      }
      qs2::qs_read(file.path(dir_data, "memory.qs2"))
    } else {
      readRDS(file.path(dir_data, "memory.rds"))
    }

    S7::prop(object, "dims") <- as.integer(rust_con$get_shape())
    # objects saved before a slot was added come back without it
    S7::prop(object, "sc_map") <- .migrate_sc_map(saved_data$sc_map)
    S7::prop(object, "sc_cache") <- .migrate_sc_cache(saved_data$sc_cache)
    S7::prop(object, "adt_cache") <- saved_data$adt_cache
    S7::prop(object, "atac_cache") <- saved_data$atac_cache
    S7::prop(object, "adt_counts") <- saved_data$adt_counts
    S7::prop(object, "other_data") <- saved_data$other_data

    if (
      (length(get_cell_names(object)) != rust_con$get_shape()[1]) |
        (length(get_gene_names(object)) != rust_con$get_shape()[2])
    ) {
      stop(paste(
        "The dimensions of the found data do not agree with the Rust data."
      ))
    }

    # ADT is keyed by barcode and frozen at add time, so it can hold more
    # cells than survive later QC; only kept cells missing from it matter
    adt <- S7::prop(object, "adt_counts")
    if (!is.null(adt)) {
      n_missing <- sum(
        !get_cell_names(object, filtered = TRUE) %in% rownames(adt$raw_counts)
      )
      if (n_missing > 0L) {
        warning(sprintf(
          "%i cells to keep are not present in the stored ADT counts.",
          n_missing
        ))
      }
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
