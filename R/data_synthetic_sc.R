# single cell synthetic data ---------------------------------------------------

## data generation -------------------------------------------------------------

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
  checkmate::qassert(seed, "I1")
  assertScSyntheticData(syn_data_params)

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
      cell_configs = marker_genes,
      batch_effect_strength = batch_effect_strength,
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

## data saving -----------------------------------------------------------------

### write h5ad type formats ----------------------------------------------------

#' Helper function to write data to h5ad format
#'
#' @description This is a helper to write synthetic data to h5ad file.
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

#' Download PBMC3K data from 10x Genomics
#'
#' @description
#' This function downloads the PBMC3K dataset from 10x Genomics and extracts
#' it to a temporary directory. It returns the path to the extracted data.
#'
#' @returns String. The path to the extracted PBMC3K data.
#' @export
#'
download_pbmc3k <- function() {
  temp_dir <- tempdir()
  dest_file <- file.path(temp_dir, "pbmc3k.tar.gz")
  url <- "https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"

  download.file(url, dest_file, mode = "wb", quiet = TRUE)
  untar(dest_file, exdir = temp_dir)

  file.path(temp_dir, "filtered_gene_bc_matrices", "hg19")
}
