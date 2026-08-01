# helpers ----------------------------------------------------------------------

## file discovery --------------------------------------------------------------

#' Locate exactly one file in a directory
#'
#' @description
#' Case-insensitive [list.files()] wrapper that insists on exactly one hit and
#' errors otherwise. Loaders use this to resolve the 10x trio and the Visium
#' `spatial/` files without hard-coding names that shift between Cell Ranger
#' and Space Ranger versions.
#'
#' @param dir String. Directory to search in.
#' @param pat String. Regular expression matched against the file names.
#'
#' @return String. Full path to the single matching file.
#'
#' @keywords internal
locate <- function(dir, pat) {
  # checks
  checkmate::assertDirectoryExists(dir)
  checkmate::qassert(pat, "S1")

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

  return(files)
}

#' Decompress a gzipped file into a temporary directory
#'
#' @description
#' Streams a `.gz` file through 8 MB chunks into `temp_dir`. The Rust mtx
#' readers want a plain file on disk, and the 10x outputs ship compressed. The
#' caller owns the result and should `unlink()` it when done.
#'
#' @param path String. Path to the gzipped file.
#' @param temp_dir String. Existing directory the decompressed file lands in.
#'
#' @return String. Path to the decompressed file.
#'
#' @keywords internal
gunzip_to_temp <- function(path, temp_dir) {
  # checks
  checkmate::assertFileExists(path)
  checkmate::assertDirectoryExists(temp_dir)

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

  return(out_path)
}

## conversion ------------------------------------------------------------------

#' Convert legacy v2 single-cell data files to v3 format
#'
#' @description Takes a directory containing legacy \code{counts_cells.bin} and
#'   \code{counts_genes.bin} files and converts them to the v3 format in place.
#'
#' @param dir_path Character scalar. Path to the directory containing the
#'   legacy binary files.
#' @param verbose Logical scalar. Whether to print progress messages.
#'   Default \code{TRUE}.
#'
#' @return Invisible \code{NULL}. Called for its side effects.
#'
#' @export
sc_old_file_conversion <- function(dir_path, verbose = TRUE) {
  checkmate::assert_string(dir_path)
  checkmate::assert_directory_exists(dir_path)
  checkmate::assert_flag(verbose)

  cell_file <- file.path(dir_path, "counts_cells.bin")
  gene_file <- file.path(dir_path, "counts_genes.bin")
  checkmate::assert_file_exists(cell_file)
  checkmate::assert_file_exists(gene_file)

  cell_tmp <- file.path(dir_path, "counts_cells_tmp.bin")
  gene_tmp <- file.path(dir_path, "counts_genes_tmp.bin")

  rs_data_v2_3_conversion(
    cell_input = cell_file,
    cell_output = cell_tmp,
    gene_input = gene_file,
    gene_output = gene_tmp,
    verbose = verbose
  )

  file.remove(cell_file)
  file.remove(gene_file)
  file.rename(cell_tmp, cell_file)
  file.rename(gene_tmp, gene_file)

  if (verbose) {
    message("Conversion complete.")
  }

  invisible(NULL)
}
