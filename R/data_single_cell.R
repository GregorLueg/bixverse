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

### dialogue -------------------------------------------------------------------

#' Single cell test data with a planted multicellular programme
#'
#' @description
#' Generates synthetic data DIALOGUE should be able to solve. Every cell type
#' gets its own noise and its own sample-level nuisance factors; only the first
#' feature column and the planted genes carry a shared per-sample latent, so
#' anything found beyond that is spurious.
#'
#' @details
#' Same generator the Rust integration tests use, so the planted structure
#' cannot drift between the two suites. What differs is the scale: R takes the
#' raw count layer and lets the normal ingestion path log-normalise it, where
#' the Rust tests feed the planted layer straight in.
#'
#' Cells are laid out contiguously by cell type and, within a cell type, by
#' sample, so every cell type sees every sample. That is the easy case for the
#' method and is what a fixture should be.
#'
#' @param syn_data_params List. Contains the parameters for the generation of
#' the synthetic data, see: [bixverse::params_sc_synthetic_dialogue()].
#' @param seed Integer. The seed for the generation of the synthetic data.
#'
#' @returns List with the following items
#' \itemize{
#'   \item counts - dgRMatrix with cells x genes of raw counts.
#'   \item obs - data.table with `cell_id`, `cell_grp`, `sample_id` and
#'   `cell_quality`. The last is pure noise, independent of the planted
#'   programme: hand it to `dialogue_sc(quality_col = "cell_quality")` so the
#'   covariate does not carry the signal you are trying to find.
#'   \item var - data.table with the gene information.
#'   \item features - Named list of numeric matrices, one per cell type, with
#'   the cell identifiers as row names. Hand this straight to
#'   [bixverse::dialogue_sc()].
#'   \item latent - Numeric vector. The per-sample latent the planted
#'   programme follows, named by sample.
#'   \item planted - Named list of character vectors. The planted gene
#'   identifiers per cell type.
#' }
#'
#' @references Jerby-Arnon & Regev, Nature Biotechnology, 2022
#'
#' @export
#'
#' @keywords internal
generate_dialogue_test_data <- function(
  syn_data_params = params_sc_synthetic_dialogue(),
  seed = 42L
) {
  # checks
  assertScSyntheticDialogue(syn_data_params)
  checkmate::qassert(seed, "I1")

  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop(
      "Package 'Matrix' needed for this function to work. Please install it.",
      call. = FALSE
    )
  }

  data <- with(
    syn_data_params,
    rs_synthetic_sc_dialogue_data(
      n_samples = n_samples,
      cells_per_sample = cells_per_sample,
      n_cell_types = n_cell_types,
      n_features = n_features,
      n_sample_features = n_sample_features,
      n_genes = n_genes,
      n_planted = n_planted,
      seed = seed
    )
  )

  n_cells <- data$nrow
  n_genes <- data$ncol

  cell_ids <- sprintf(
    sprintf("cell_%%0%dd", nchar(as.character(n_cells))),
    seq_len(n_cells)
  )
  gene_ids <- sprintf(
    sprintf("gene_%%0%dd", nchar(as.character(n_genes))),
    seq_len(n_genes)
  )
  cell_types <- sprintf("cell_type_%i", seq_along(data$cell_type_indices))
  sample_ids <- sprintf("sample_%02d", data$sample_ids + 1L)

  counts <- new(
    "dgRMatrix",
    p = as.integer(data$indptr),
    x = as.numeric(data$data),
    j = as.integer(data$indices),
    Dim = as.integer(c(n_cells, n_genes)),
    Dimnames = list(cell_ids, gene_ids)
  )

  cell_grp <- character(n_cells)
  for (i in seq_along(data$cell_type_indices)) {
    cell_grp[data$cell_type_indices[[i]] + 1L] <- cell_types[i]
  }

  obs <- data.table::data.table(
    cell_id = cell_ids,
    cell_grp = cell_grp,
    sample_id = sample_ids,
    cell_quality = data$quality
  )

  var <- data.table::data.table(
    gene_id = gene_ids,
    ensembl_id = sprintf(
      sprintf("ens_%%0%dd", nchar(as.character(n_genes))),
      seq_len(n_genes)
    )
  )

  features <- purrr::imap(data$features, \(mat, i) {
    rownames(mat) <- cell_ids[data$cell_type_indices[[i]] + 1L]
    colnames(mat) <- sprintf("feature_%02d", seq_len(ncol(mat)))
    mat
  })
  names(features) <- cell_types

  planted <- purrr::map(data$planted, \(idx) gene_ids[idx + 1L])
  names(planted) <- cell_types

  latent <- data$latent
  names(latent) <- sprintf("sample_%02d", seq_along(latent))

  list(
    counts = counts,
    obs = obs,
    var = var,
    features = features,
    latent = latent,
    planted = planted
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

#' Resolve the download URLs for an example data file
#'
#' @description
#' GitHub releases first, Zenodo second. Zenodo holds the archival copy with
#' the DOI, but it is slow and drops connections often enough to fail a
#' vignette build outright.
#'
#' @param file String. File name, identical in both locations.
#' @param zenodo_record Optional string. Zenodo record holding the fallback
#' copy. `NULL` for datasets that only live in the GitHub release, which is the
#' case until they are archived.
#'
#' @returns Character vector of URLs, in the order they should be tried.
#'
#' @keywords internal
.data_urls <- function(file, zenodo_record = NULL) {
  checkmate::assertString(file)
  checkmate::assertString(zenodo_record, null.ok = TRUE)

  gh <- sprintf(
    "https://github.com/GregorLueg/bixverse-data/releases/download/%s/%s",
    .BIXVERSE_DATA_TAG,
    file
  )

  # not everything has an archival copy yet, and a dataset that only lives in
  # the release is still perfectly downloadable
  if (is.null(zenodo_record)) {
    return(gh)
  }

  c(
    gh,
    sprintf(
      "https://zenodo.org/records/%s/files/%s?download=1",
      zenodo_record,
      file
    )
  )
}

#' Release tag of the bixverse-data repository
#'
#' @description
#' Pinned rather than moving, so a given bixverse version always resolves the
#' same data.
#'
#' @keywords internal
.BIXVERSE_DATA_TAG <- "v1.0.0"

#' Download a file, retrying on failure
#'
#' @description
#' Zenodo drops connections often enough that a single `download.file()` call
#' is not reliable, which shows up as `cannot open URL` mid-vignette. Retries
#' with a linear backoff.
#'
#' @param url String. The URL to fetch.
#' @param dest_file String. Path to write the file to.
#' @param quiet Boolean. If the download shall be quiet.
#' @param tries Integer. Attempts before giving up.
#'
#' @returns Invisibly, `dest_file`.
#'
#' @keywords internal
.download_with_retry <- function(urls, dest_file, quiet = FALSE, tries = 4L) {
  checkmate::assertCharacter(urls, min.len = 1, any.missing = FALSE)
  checkmate::assertString(dest_file)
  checkmate::assertFlag(quiet)
  checkmate::assertCount(tries, positive = TRUE)

  old_timeout <- getOption("timeout")
  options(timeout = max(300, old_timeout))
  on.exit(options(timeout = old_timeout))

  for (i in seq_len(tries)) {
    for (url in urls) {
      ok <- tryCatch(
        {
          download.file(url, dest_file, mode = "wb", quiet = quiet)
          TRUE
        },
        error = function(e) FALSE,
        warning = function(w) FALSE
      )

      if (ok && file.exists(dest_file) && file.size(dest_file) > 0) {
        return(invisible(dest_file))
      }

      unlink(dest_file)
    }

    if (i < tries) {
      cli::cli_alert_warning("Download failed ({i}/{tries}), retrying.")
      Sys.sleep(5 * i)
    }
  }

  cli::cli_abort("Download failed after {tries} attempts: {.url {urls[[1]]}}")
}


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
  urls <- .data_urls("pbmc3k_filtered_gene_bc_matrices.tar.gz", "20977604")

  .download_with_retry(urls, dest_file, quiet = quiet)
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
  urls <- .data_urls("pmbc-8k.tar.gz", "20977604")

  .download_with_retry(urls, dest_file, quiet = quiet)
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
  urls <- .data_urls("demuxlet_PBMCs.tar.gz", "20977604")

  .download_with_retry(urls, dest_file, quiet = quiet)
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
  urls <- .data_urls("pbmc_batches.tar.gz", "20977604")

  .download_with_retry(urls, dest_file, quiet = quiet)
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
  urls <- .data_urls("cd34_multiome_rna.h5ad.gz", "20977604")

  .download_with_retry(urls, dest_file, quiet = quiet)
  R.utils::gunzip(dest_file, remove = TRUE)

  file.path(temp_dir, "cd34_multiome_rna.h5ad")
}

### dialogue ulcerative colitis ------------------------------------------------

#' Download the ulcerative colitis example data for DIALOGUE
#'
#' @description
#' This function downloads the ulcerative colitis subset used by the DIALOGUE
#' paper into the temporary directory. 5,374 cells across five cell subtypes
#' and 30 donors, with the inflammation status per biopsy. Enough samples per
#' cell type for [bixverse::dialogue_sc()], which most single sample example
#' data sets are not.
#'
#' @details
#' The published matrix holds `log2(TPM / 10 + 1)`, which is what Smillie, et
#' al. released; the raw counts were never published. The file served here has
#' been linearised back to the TPM/10 scale and rounded, so it loads through
#' [bixverse::load_h5ad()] like any count matrix and gets bixverse's own log
#' CPM on the way in. See `data-raw/dialogue_uc.R` for the preparation.
#'
#' @param quiet Boolean. If the download shall be quiet.
#'
#' @returns String. The path to the ulcerative colitis h5ad file.
#'
#' @export
#'
#' @references Smillie, et al., Cell, 2019; Jerby-Arnon and Regev, Nat.
#' Biotechnol., 2022
download_dialogue_uc <- function(quiet = FALSE) {
  temp_dir <- tempdir()
  dest_file <- file.path(temp_dir, "dialogue_uc.h5ad.gz")
  urls <- .data_urls("dialogue_uc.h5ad.gz", "22105703")

  .download_with_retry(urls, dest_file, quiet = quiet)
  R.utils::gunzip(dest_file, remove = TRUE)

  file.path(temp_dir, "dialogue_uc.h5ad")
}

### marrow cd34 example data set -----------------------------------------------

#' Download the marrow CD34 example data from Palantir
#'
#' @description
#' This function downloads the bone marrow CD34 data set from the Palantir
#' paper into the temporary directory.
#'
#' @param quiet Boolean. If the download shall be quiet.
#'
#' @returns String. The path to the marrow CD34 data set.
#'
#' @export
#'
#' @references Setty, et al., Nat. Biotechnol., 2019
download_marrow_cd34 <- function(quiet = FALSE) {
  temp_dir <- tempdir()
  dest_file <- file.path(temp_dir, "marrow_sample_scseq_counts.h5ad.gz")
  urls <- .data_urls("marrow_sample_scseq_counts.h5ad.gz", "21892320")

  .download_with_retry(urls, dest_file, quiet = quiet)
  R.utils::gunzip(dest_file, remove = TRUE)

  file.path(temp_dir, "marrow_sample_scseq_counts.h5ad")
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
  urls <- .data_urls("10k_Human_PBMC_TotalSeqB.h5", "20977604")

  .download_with_retry(urls, dest_file, quiet = quiet)

  file.path(temp_dir, "10k_Human_PBMC_TotalSeqB.h5")
}

### kang pbmc ------------------------------------------------------------------

#' Download the Kang, et al. IFN-beta stimulated PBMC data
#'
#' @description
#' Downloads the Kang, et al. demultiplexed PBMC experiment as a
#' `SingleCellExperiment`. 8 lupus donors, PBMCs, control against IFN-beta
#' stimulated in vitro, roughly 29k cells before quality control.
#'
#' The design is paired within donor, which is what makes it the example for
#' [bixverse::nebula_sc()]: the donor is a genuine random effect and every donor
#' contributes to both arms. Load it with [bixverse::load_sce()].
#'
#' @details
#' Sourced from the `muscData` Bioconductor package, itself from GEO accession
#' `GSE96583`. See `data-raw/zenodo_sce_datasets.R` for the preparation. The
#' `multiplets` column marks doublets and ambiguous droplets, which are worth
#' dropping before any modelling.
#'
#' Served from the `bixverse-data` GitHub release. There is no Zenodo fallback
#' for this one yet.
#'
#' @param quiet Boolean. If the download shall be quiet.
#'
#' @returns String. The path to the qs2 file holding the
#' `SingleCellExperiment`.
#'
#' @export
#'
#' @references Kang, et al., Nat. Biotechnol., 2018
download_kang_pbmc <- function(quiet = FALSE) {
  temp_dir <- tempdir()
  dest_file <- file.path(temp_dir, "kang18_8vs8.qs2")
  urls <- .data_urls("kang18_8vs8.qs2")

  .download_with_retry(urls, dest_file, quiet = quiet)

  dest_file
}

### ageing thymus --------------------------------------------------------------

#' Download the Baran-Gale, et al. ageing thymus data
#'
#' @description
#' Downloads the mouse ageing thymus droplet experiment as a
#' `SingleCellExperiment`. Roughly 69k thymic epithelial cells across three
#' ages, with real change in cell type proportions between them.
#'
#' That compositional change is what makes it the example for
#' [bixverse::get_miloR_abundances_sc()] and [bixverse::meld_sc()]. A
#' stimulation experiment moves cells in embedding space without moving the
#' proportions, so almost every neighbourhood comes out significant and the
#' result says nothing. Load it with [bixverse::load_sce()].
#'
#' @details
#' Sourced from the `MouseThymusAgeing` Bioconductor package. See
#' `data-raw/zenodo_sce_datasets.R` for the preparation. The object carries no
#' gene names, so [bixverse::load_sce()] falls back to the Ensembl identifiers
#' in the `rowData`, and 336 genes with no annotation at all get a generated
#' identifier.
#'
#' Served from the `bixverse-data` GitHub release. There is no Zenodo fallback
#' for this one yet. At 430 MB it is the largest of the example datasets, so
#' expect the first call to take a while.
#'
#' @param quiet Boolean. If the download shall be quiet.
#'
#' @returns String. The path to the qs2 file holding the
#' `SingleCellExperiment`.
#'
#' @export
#'
#' @references Baran-Gale, et al., Development, 2020
download_thymus_ageing <- function(quiet = FALSE) {
  temp_dir <- tempdir()
  dest_file <- file.path(temp_dir, "thymus_ageing_droplet.qs2")
  urls <- .data_urls("thymus_ageing_droplet.qs2")

  .download_with_retry(urls, dest_file, quiet = quiet)

  dest_file
}
