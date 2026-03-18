# helpers ----------------------------------------------------------------------

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
