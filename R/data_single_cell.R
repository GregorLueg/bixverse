# single cell synthetic data ---------------------------------------------------

## data generation -------------------------------------------------------------

### rna ------------------------------------------------------------------------

#' Single cell test data
#'
#' @description
#' This function generates synthetic data for single cell test purposes. These
#' data can be used for testing functionality of various single cell functions.
#'
#' @param syn_data_params List. Contains the parameters for the generation of
#' synthetic data, see: [bixverse::params_sc_synthetic_data()]. Has the
#' following elements:
#' \itemize{
#'   \item n_cells - Integer. Number of cells.
#'   \item n_genes - Integer. Number of genes.
#'   \item marker_genes - List. A nested list that indicates which gene indices
#'   are markers for which cell.
#'   \item n_batches - Integer. Number of batches.
#'   \item batch_effect_strength - String. Indicates the strength of the batch
#'   effect to add.
#' }
#' @param seed Integer. The seed for the generation of the seed data.
#'
#' @returns List with the following items
#' \itemize{
#'   \item counts - dgRMatrix with cells x genes.
#'   \item obs - data.table that contains the cell information.
#'   \item var - data.table that contains the var information.
#' }
#'
#' @export
generate_single_cell_test_data <- function(
  syn_data_params = params_sc_synthetic_data(),
  seed = 42L
) {
  # checks
  assertScSyntheticData(syn_data_params)
  checkmate::qassert(seed, "I1")

  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop(
      "Package 'Matrix' needed for this function to work. Please install it.",
      call. = FALSE
    )
  }

  data <- with(
    syn_data_params,
    rs_synthetic_sc_data_with_cell_types(
      n_cells = n_cells,
      n_genes = n_genes,
      n_batches = n_batches,
      n_samples = n_samples,
      cell_configs = marker_genes,
      batch_effect_strength = batch_effect_strength,
      sample_bias = sample_bias,
      seed = seed
    )
  )

  counts <- with(
    syn_data_params,
    new(
      "dgRMatrix",
      p = as.integer(data$indptr),
      x = as.numeric(data$data),
      j = as.integer(data$indices),
      Dim = as.integer(c(n_cells, n_genes))
    )
  )

  n_digits <- nchar(as.character(syn_data_params$n_cells))
  format_str <- sprintf("cell_%%0%dd", n_digits)
  rownames(counts) <- sprintf(
    format_str,
    1:syn_data_params$n_cells
  )

  obs <- data.table(
    cell_id = sprintf(
      format_str,
      1:syn_data_params$n_cells
    ),
    cell_grp = sprintf("cell_type_%i", data$cell_type_indices + 1),
    batch_index = data$batch_indices + 1
  )

  if (!is.null(data$sample_indices)) {
    obs[, sample_id := sprintf("sample_%i", data$sample_indices + 1)]
  }

  n_digits <- nchar(as.character(syn_data_params$n_genes))
  format_str <- sprintf("gene_%%0%dd", n_digits)
  format_str_2 <- sprintf("ens_%%0%dd", n_digits)
  colnames(counts) <- sprintf(
    format_str,
    1:syn_data_params$n_genes
  )

  var <- data.table(
    gene_id = sprintf(
      format_str,
      1:syn_data_params$n_genes
    ),
    ensembl_id = sprintf(
      format_str_2,
      1:syn_data_params$n_genes
    )
  )

  res <- list(
    counts = counts,
    obs = obs,
    var = var
  )
}

### adt ------------------------------------------------------------------------

#' Single cell test data (ADT)
#'
#' @description
#' This function generates synthetic ADT counts for single cell test purposes.
#' These data can be used for testing multi-modal functionality of various
#' single cell functions. Pairs cell-for-cell with
#' [bixverse::generate_single_cell_test_data()] when generated with matching
#' `n_cells`, number of cell types and `n_batches`.
#'
#' @param syn_data_params List. Contains the parameters for the generation of
#' synthetic ADT data, see: [bixverse::params_sc_synthetic_data_adt()]. Has the
#' following elements:
#' \itemize{
#'   \item n_cells - Integer. Number of cells.
#'   \item n_proteins - Integer. Number of proteins.
#'   \item marker_genes - List. A nested list that indicates which protein
#'   indices are markers for which cell type.
#'   \item n_batches - Integer. Number of batches.
#'   \item isotype_controls - Integer vector. Column indices (0-based) of the
#'   isotype controls.
#'   \item batch_effect_strength - String. Indicates the strength of the batch
#'   effect to add.
#' }
#' @param seed Integer. The seed for the generation of the synthetic data.
#'
#' @returns List with the following items
#' \itemize{
#'   \item counts - Numeric matrix with cells x proteins.
#'   \item obs - data.table that contains the cell information.
#'   \item var - data.table that contains the protein information.
#' }
#'
#' @export
#'
#' @keywords internal
generate_single_cell_test_data_adt <- function(
  syn_data_params = params_sc_synthetic_data_adt(),
  seed = 42L
) {
  # checks
  assertScSyntheticDataAdt(syn_data_params)
  checkmate::qassert(seed, "I1")

  data <- with(
    syn_data_params,
    rs_synthetic_sc_adt_with_cell_types(
      n_cells = n_cells,
      n_proteins = n_proteins,
      n_batches = n_batches,
      isotype_controls = isotype_controls,
      cell_configs = marker_genes,
      batch_effect_strength = batch_effect_strength,
      seed = seed
    )
  )

  # the Rust side returns a flat row-major (cell-major) vector
  counts <- matrix(
    as.numeric(data$data),
    nrow = syn_data_params$n_cells,
    ncol = syn_data_params$n_proteins,
    byrow = TRUE
  )

  n_digits <- nchar(as.character(syn_data_params$n_cells))
  format_str <- sprintf("cell_%%0%dd", n_digits)
  rownames(counts) <- sprintf(format_str, 1:syn_data_params$n_cells)

  obs <- data.table(
    cell_id = sprintf(format_str, 1:syn_data_params$n_cells),
    cell_grp = sprintf("cell_type_%i", data$cell_type_indices + 1),
    batch_index = data$batch_indices + 1
  )

  n_digits <- nchar(as.character(syn_data_params$n_proteins))
  format_str <- sprintf("protein_%%0%dd", n_digits)
  colnames(counts) <- sprintf(format_str, 1:syn_data_params$n_proteins)

  isotype_ids <- sprintf(format_str, syn_data_params$isotype_controls + 1)

  var <- data.table(
    protein_id = sprintf(format_str, 1:syn_data_params$n_proteins),
    is_isotype = sprintf(format_str, 1:syn_data_params$n_proteins) %in%
      isotype_ids
  )

  res <- list(
    counts = counts,
    obs = obs,
    var = var
  )

  res
}

## data saving -----------------------------------------------------------------

### write h5ad type formats ----------------------------------------------------

#### sparse --------------------------------------------------------------------

#' Helper function to write data to h5ad format
#'
#' @description This is a helper to write synthetic data to h5ad file. This
#' version will write the data into the common compressed sparse data format.
#'
#' @param f_path String. The filepath to which to save the data
#' @param counts Sparse matrix. Needs to be of class `dgRMatrix` or
#' `dgCMatrix`.
#' @param obs data.table. The observations. Needs to have
#' `nrow(obs) == nrow(counts)`.
#' @param var data.table. The variable data. Needs to have
#' `ncol(var) == ncol(counts)`.
#' @param overwrite Boolean. Shall any found h5ad file be overwritten.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @return Returns invisible
#'
#' @export
write_h5ad_sc <- function(
  f_path,
  counts,
  obs,
  var,
  overwrite = TRUE,
  .verbose = TRUE
) {
  # checks
  checkmate::assertPathForOutput(f_path, overwrite = TRUE)
  checkmate::assert(
    checkmate::testClass(counts, "dgRMatrix"),
    checkmate::testClass(counts, "dgCMatrix")
  )
  checkmate::assertDataTable(
    obs,
    min.rows = nrow(counts),
    max.rows = nrow(counts)
  )
  checkmate::assertDataTable(
    var,
    min.rows = ncol(counts),
    max.rows = ncol(counts)
  )
  checkmate::qassert(overwrite, "B1")
  checkmate::qassert(.verbose, "B1")

  if (file.exists(f_path) & !overwrite) {
    stop("The h5ad file already exists and overwrite = FALSE.")
  } else if (file.exists(f_path)) {
    file.remove(f_path)
  }

  rhdf5::h5createFile(f_path)

  if (.verbose) {
    message("Writing the counts to h5ad.")
  }
  # Write X (sparse matrix)
  rhdf5::h5createGroup(f_path, "X")

  if (inherits(counts, "dgCMatrix")) {
    # CSC format
    rhdf5::h5write(counts@x, f_path, "X/data")
    rhdf5::h5write(counts@i, f_path, "X/indices")
    rhdf5::h5write(counts@p, f_path, "X/indptr")
  } else if (inherits(counts, "dgRMatrix")) {
    # CSR format
    rhdf5::h5write(counts@x, f_path, "X/data")
    rhdf5::h5write(counts@j, f_path, "X/indices")
    rhdf5::h5write(counts@p, f_path, "X/indptr")
  }

  if (.verbose) {
    message("Writing the obs to h5ad.")
  }

  rhdf5::h5createGroup(f_path, "obs")
  rhdf5::h5write(obs[[1]], f_path, "obs/_index")
  for (col in names(obs)[-1]) {
    rhdf5::h5write(obs[[col]], f_path, paste0("obs/", col))
  }

  if (.verbose) {
    message("Writing the var to h5ad.")
  }

  rhdf5::h5createGroup(f_path, "var")
  rhdf5::h5write(var[[1]], f_path, "var/_index")
  for (col in names(var)[-1]) {
    rhdf5::h5write(var[[col]], f_path, paste0("var/", col))
  }

  rhdf5::h5closeAll()

  invisible()
}

#### dense ---------------------------------------------------------------------

#' Helper function to write data to a dense h5ad file
#'
#' @description Helper to write synthetic data to h5ad with `/X` stored as a
#' dense 2D dataset (cells x genes). Useful for exercising the dense ingestion
#' path.
#'
#' @param f_path String. The filepath to which to save the data
#' @param counts Matrix or sparse matrix; sparse input is densified.
#' @param obs data.table. The observations. Needs `nrow(obs) == nrow(counts)`.
#' @param var data.table. The variable data. Needs `nrow(var) == ncol(counts)`.
#' @param overwrite Boolean. Shall any found h5ad file be overwritten.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @return Returns invisible
#'
#' @export
write_h5ad_sc_dense <- function(
  f_path,
  counts,
  obs,
  var,
  overwrite = TRUE,
  .verbose = TRUE
) {
  checkmate::assertPathForOutput(f_path, overwrite = TRUE)
  checkmate::assertDataTable(
    obs,
    min.rows = nrow(counts),
    max.rows = nrow(counts)
  )
  checkmate::assertDataTable(
    var,
    min.rows = ncol(counts),
    max.rows = ncol(counts)
  )
  checkmate::qassert(overwrite, "B1")
  checkmate::qassert(.verbose, "B1")

  if (file.exists(f_path) && !overwrite) {
    stop("The h5ad file already exists and overwrite = FALSE.")
  } else if (file.exists(f_path)) {
    file.remove(f_path)
  }

  rhdf5::h5createFile(f_path)

  counts_dense <- as.matrix(counts)

  if (.verbose) {
    message("Writing the dense counts to h5ad.")
  }
  # Dense X as a 2D dataset directly under /X
  rhdf5::h5write(counts_dense, f_path, "X", native = TRUE)

  if (.verbose) {
    message("Writing the obs to h5ad.")
  }
  rhdf5::h5createGroup(f_path, "obs")
  rhdf5::h5write(obs[[1]], f_path, "obs/_index")
  for (col in names(obs)[-1]) {
    rhdf5::h5write(obs[[col]], f_path, paste0("obs/", col))
  }

  if (.verbose) {
    message("Writing the var to h5ad.")
  }
  rhdf5::h5createGroup(f_path, "var")
  rhdf5::h5write(var[[1]], f_path, "var/_index")
  for (col in names(var)[-1]) {
    rhdf5::h5write(var[[col]], f_path, paste0("var/", col))
  }

  rhdf5::h5closeAll()

  invisible()
}

### write cell ranger output ---------------------------------------------------

#' Helper function to write data to a cell ranger like output
#'
#' @description This is a helper to write synthetic data to cell ranger like
#' output, i.e., an .mtx file, an barcodes.csv (or .tsv) and a features.csv (or
#' .tsv).
#'
#' @param f_path String. The filepath to which to save the data
#' @param counts Sparse matrix. Needs to be of class `dgRMatrix` or
#' `dgCMatrix`.
#' @param obs data.table. The observations. Needs to have
#' `nrow(obs) == nrow(counts)`.
#' @param var data.table. The variable data. Needs to have
#' `ncol(var) == ncol(counts)`.
#' @param format_type String. One of `c("csv", "tsv")`. Shall the data be
#' saved in TSV or CSV.
#' @param rows String. One of `c("cells", "genes")`. Shall the rows represent
#' cells or genes in the .mtx file.
#' @param overwrite Boolean. Shall any found h5ad file be overwritten.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @return Returns invisible
#'
#' @export
write_cellranger_output <- function(
  f_path,
  counts,
  obs,
  var,
  format_type = c("csv", "tsv"),
  rows = c("cells", "genes"),
  overwrite = TRUE,
  .verbose = TRUE
) {
  format_type <- match.arg(format_type)
  rows <- match.arg(rows)

  # checks
  checkmate::assertPathForOutput(f_path, overwrite = TRUE)
  checkmate::assert(
    checkmate::testClass(counts, "dgRMatrix"),
    checkmate::testClass(counts, "dgCMatrix")
  )
  checkmate::assertDataTable(
    obs,
    min.rows = nrow(counts),
    max.rows = nrow(counts)
  )
  checkmate::assertDataTable(
    var,
    min.rows = ncol(counts),
    max.rows = ncol(counts)
  )
  checkmate::assertChoice(format_type, c("csv", "tsv"))
  checkmate::assertChoice(rows, c("cells", "genes"))

  f_path_mtx <- file.path(f_path, "mat.mtx")
  f_path_obs <- file.path(f_path, sprintf("barcodes.%s", format_type))
  f_path_var <- file.path(f_path, sprintf("features.%s", format_type))

  if (
    (file.exists(f_path_mtx) ||
      file.exists(f_path_obs) ||
      file.exists(f_path_var)) &
      !overwrite
  ) {
    stop("The to be written files already exist and overwrite = FALSE.")
  }

  # extract triplets from dgRMatrix
  row_idx <- rep(seq_len(nrow(counts)), diff(counts@p))
  col_idx <- counts@j + 1L
  values <- counts@x

  if (rows == "cells") {
    rows_dim <- nrow(counts)
    cols_dim <- ncol(counts)
  } else {
    # swap for genes x cells
    rows_dim <- ncol(counts)
    cols_dim <- nrow(counts)
    temp <- row_idx
    row_idx <- col_idx
    col_idx <- temp
  }

  if (.verbose) {
    message("Writing the .mtx file to disk.")
  }

  # write matrix.mtx
  mtx_file <- file.path(f_path_mtx)
  con <- file(mtx_file, "w")
  writeLines("%%MatrixMarket matrix coordinate integer general", con)
  writeLines(
    sprintf(
      "%%Rows=%s (%d), Cols=%s (%d)",
      ifelse(rows == "cells", "cells", "genes"),
      rows_dim,
      ifelse(rows == "cells", "genes", "cells"),
      cols_dim
    ),
    con
  )
  writeLines(sprintf("%d %d %d", rows_dim, cols_dim, length(values)), con)
  writeLines(sprintf("%d %d %d", row_idx, col_idx, as.integer(values)), con)
  close(con)

  if (.verbose) {
    message("Writing the obs and var files to disk.")
  }

  # write barcodes / obs
  data.table::fwrite(
    x = obs,
    file = f_path_obs,
    sep = ifelse(format_type == "csv", ",", "\t")
  )

  # write features / var
  data.table::fwrite(
    x = var,
    file = f_path_var,
    sep = ifelse(format_type == "csv", ",", "\t")
  )

  invisible()
}

### write 10x h5 files ---------------------------------------------------------

#' Helper function to write data to a 10x CellRanger-style h5 file
#'
#' @param f_path String. Output path.
#' @param counts Sparse matrix (`dgRMatrix` or `dgCMatrix`), cells x features.
#' @param barcodes Character. Cell barcodes, length `nrow(counts)`.
#' @param features data.table with `id` and `name` of length `ncol(counts)`. For
#' v3 may include `feature_type`; defaults to `"Gene Expression"` if absent.
#' @param version One of `"v3"` or `"v2"`.
#' @param overwrite Boolean.
#'
#' @return Invisible.
#'
#' @export
#'
#' @keywords internal
write_tenx_h5_sc <- function(
  f_path,
  counts,
  barcodes,
  features,
  version = c("v3", "v2"),
  overwrite = TRUE
) {
  version <- match.arg(version)

  checkmate::assertPathForOutput(f_path, overwrite = TRUE)
  checkmate::assert(
    checkmate::testClass(counts, "dgRMatrix"),
    checkmate::testClass(counts, "dgCMatrix")
  )
  checkmate::assertCharacter(barcodes, len = nrow(counts))
  checkmate::assertDataTable(features, nrows = ncol(counts))
  checkmate::assertSubset(c("id", "name"), names(features))
  checkmate::qassert(overwrite, "B1")

  if (file.exists(f_path) && !overwrite) {
    stop("The h5 file already exists and overwrite = FALSE.")
  } else if (file.exists(f_path)) {
    file.remove(f_path)
  }

  # cells x features dgRMatrix has the same slot layout as a
  # features x cells dgCMatrix (which is what 10x stores).
  counts_r <- if (inherits(counts, "dgRMatrix")) {
    counts
  } else {
    as(counts, "RsparseMatrix")
  }

  indptr <- as.integer(counts_r@p)
  indices <- as.integer(counts_r@j)
  data <- as.integer(counts_r@x)
  shape <- as.integer(c(ncol(counts), nrow(counts)))

  rhdf5::h5createFile(f_path)
  on.exit(tryCatch(rhdf5::h5closeAll(), error = function(e) invisible()))

  if (version == "v3") {
    rhdf5::h5createGroup(f_path, "matrix")
    rhdf5::h5createGroup(f_path, "matrix/features")

    rhdf5::h5write(data, f_path, "matrix/data")
    rhdf5::h5write(indices, f_path, "matrix/indices")
    rhdf5::h5write(indptr, f_path, "matrix/indptr")
    rhdf5::h5write(shape, f_path, "matrix/shape")
    rhdf5::h5write(barcodes, f_path, "matrix/barcodes")
    rhdf5::h5write(as.character(features$id), f_path, "matrix/features/id")
    rhdf5::h5write(
      as.character(features$name),
      f_path,
      "matrix/features/name"
    )

    feature_type <- if ("feature_type" %in% names(features)) {
      as.character(features$feature_type)
    } else {
      rep("Gene Expression", nrow(features))
    }
    rhdf5::h5write(feature_type, f_path, "matrix/features/feature_type")
  } else {
    rhdf5::h5write(data, f_path, "data")
    rhdf5::h5write(indices, f_path, "indices")
    rhdf5::h5write(indptr, f_path, "indptr")
    rhdf5::h5write(shape, f_path, "shape")
    rhdf5::h5write(barcodes, f_path, "barcodes")
    rhdf5::h5write(as.character(features$id), f_path, "genes")
    rhdf5::h5write(as.character(features$name), f_path, "gene_names")
  }

  invisible()
}

## example data sets -----------------------------------------------------------

### pbmc3k ---------------------------------------------------------------------

#' Download PBMC3K data from Zenodo
#'
#' @description
#' This function downloads the PBMC3K dataset from 10x Genomics (uploaded on
#' Zenodo) and extracts it and returns the (temporary) paths.
#'
#' @param quiet Boolean. If the download shall be quiet.
#'
#' @returns String. The path to the extracted PBMC3K data.
#'
#' @export
download_pbmc3k <- function(quiet = FALSE) {
  temp_dir <- tempdir()
  dest_file <- file.path(temp_dir, "pbmc3k.tar.gz")
  url <- "https://zenodo.org/records/20977604/files/pbmc3k_filtered_gene_bc_matrices.tar.gz?download=1"

  download.file(url, dest_file, mode = "wb", quiet = quiet)
  untar(dest_file, exdir = temp_dir)

  # add headers to genes.tsv
  data_path <- file.path(temp_dir, "filtered_gene_bc_matrices", "hg19")

  data_path
}

### pbmc8k ---------------------------------------------------------------------

#' Download PBMC8K data from Zenodo
#'
#' @description
#' This function downloads the PBMC8k dataset from 10x Genomics (uploaded on
#' Zenodo) and extracts it and returns the (temporary) paths.
#'
#' @param quiet Boolean. If the download shall be quiet.
#'
#' @returns String. The path to the extracted PBMC8K data.
#'
#' @export
download_pbmc8k <- function(quiet = FALSE) {
  temp_dir <- tempdir()
  dest_file <- file.path(temp_dir, "pmbc-8k.tar.gz")
  url <- "https://zenodo.org/records/20977604/files/pmbc-8k.tar.gz?download=1"

  download.file(url, dest_file, mode = "wb", quiet = quiet)
  untar(dest_file, exdir = temp_dir)

  # add headers to genes.tsv
  data_path <- file.path(temp_dir, "pmbc-8k")

  data_path
}

### pbmc with demuxlet ---------------------------------------------------------

#' Download PBMCs with demuxlet doublet information
#'
#' @description
#' This function downloads a PBMC data set with demuxlet information to test
#' doublet detection methods.
#'
#' @param quiet Boolean. If the download shall be quiet.
#'
#' @returns String. The path to the extracted doublet detection data.
#'
#' @export
download_demuxlet_pbmc <- function(quiet = FALSE) {
  temp_dir <- tempdir()
  dest_file <- file.path(temp_dir, "demuxlet_PBMCs.tar.gz")
  url <- "https://zenodo.org/records/20977604/files/demuxlet_PBMCs.tar.gz?download=1"

  download.file(url, dest_file, mode = "wb", quiet = quiet)
  untar(dest_file, exdir = temp_dir)

  data_path <- file.path(temp_dir, "demuxlet_PBMCs")

  data_path
}

### pbmc batches ---------------------------------------------------------------

#' Download two different PBMC data sets for batch correction testing
#'
#' @description
#' This function downloads two different h5ad files for testing batch correction
#' methods into the temporary directory.
#'
#' @param quiet Boolean. If the download shall be quiet.
#'
#' @returns String. The path to the directory with the PBMC h5ad files.
#'
#' @export
download_pbmc_batches <- function(quiet = FALSE) {
  temp_dir <- tempdir()
  dest_file <- file.path(temp_dir, "pbmc_batches.tar.gz")
  url <- "https://zenodo.org/records/20977604/files/pbmc_batches.tar.gz?download=1"

  download.file(url, dest_file, mode = "wb", quiet = quiet)
  untar(dest_file, exdir = temp_dir)

  # add headers to genes.tsv
  data_path <- file.path(temp_dir, "pbmc_batches")

  data_path
}

### cd34 example data sets -----------------------------------------------------

#' Download the CD34 example data from SEACells
#'
#' @description
#' This function downloads the CD34 data set from the SEACells paper into the
#' temperorary directory.
#'
#' @param quiet Boolean. If the download shall be quiet.
#'
#' @returns String. The path to CD34 SEACells data set.
#'
#' @export
#'
#' @references Persad, et al., Nat. Biotechnol., 2023
download_cd34_data <- function(quiet = FALSE) {
  temp_dir <- tempdir()
  dest_file <- file.path(temp_dir, "cd34_multiome_rna.h5ad.gz")
  url <- "https://zenodo.org/records/20977604/files/cd34_multiome_rna.h5ad.gz?download=1"

  download.file(url, dest_file, mode = "wb", quiet = quiet)
  R.utils::gunzip(dest_file, remove = TRUE)

  file.path(temp_dir, "cd34_multiome_rna.h5ad")
}

### pbmc totalseq --------------------------------------------------------------

#' Download the PBMC TotalSeq data with ADT counts
#'
#' @description
#' This function downloads the PBMC TotalSeq data set.
#'
#' @param quiet Boolean. If the download shall be quiet.
#'
#' @returns String. The path to the TotalSeq data.
#'
#' @export
download_pbmc_totalseq_data <- function(quiet = FALSE) {
  temp_dir <- tempdir()
  dest_file <- file.path(temp_dir, "10k_Human_PBMC_TotalSeqB.h5")
  url <- "https://zenodo.org/records/20977604/files/10k_Human_PBMC_TotalSeqB.h5?download=1"

  download.file(url, dest_file, mode = "wb", quiet = quiet)

  file.path(temp_dir, "10k_Human_PBMC_TotalSeqB.h5")
}
