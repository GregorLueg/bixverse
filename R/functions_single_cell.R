# helpers ----------------------------------------------------------------------

## cell ranger outputs ---------------------------------------------------------

#' Helper to generate cell ranger input parameters
#'
#' @param dir_data String. The directory with the Cell Ranger outputs
#'
#' @return A list based on [bixverse::params_sc_mtx_io()].
#'
#' @export
get_cell_ranger_params <- function(dir_data) {
  # checks
  checkmate::assertDirectory(dir_data)
  assertFileExists(dir_data, c("barcodes.tsv", "genes.tsv", "matrix.mtx"))

  res <- params_sc_mtx_io(
    path_mtx = path.expand(file.path(dir_data, "matrix.mtx")),
    path_obs = path.expand(file.path(dir_data, "barcodes.tsv")),
    path_var = path.expand(file.path(dir_data, "genes.tsv")),
    cells_as_rows = FALSE,
    has_hdr = FALSE
  )

  return(res)
}

## seurat assay to list --------------------------------------------------------

#' Transform Seurat raw counts into a List
#'
#' @param seurat_obj `Seurat` class. The class from which to extract the counts
#' from and transform to a list.
#'
#' @return A list with the following elements
#' \itemize{
#'   \item indptr - Index pointers of the sparse data.
#'   \item indices - Indices of the data.
#'   \item data - The underlying data.
#'   \item format - String that defines if the data is CSR or CSC.
#'   \item nrow - The number of rows.
#'   \item ncol - The number of columns.
#' }
get_seurat_counts_to_list <- function(seurat_obj) {
  # checks
  checkmate::assertClass(seurat_obj, "Seurat")

  # get the counts
  raw_counts <- seurat_obj@assays$RNA@layers$counts

  res <- list(
    indptr = raw_counts@p,
    indices = raw_counts@i,
    data = raw_counts@x,
    format = "csr",
    nrow = raw_counts@Dim[2],
    ncol = raw_counts@Dim[1]
  )

  return(res)
}
