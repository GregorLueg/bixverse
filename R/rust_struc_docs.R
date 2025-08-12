#' Single Cell Count Data Handler
#'
#' A class for efficiently reading and writing single cell count matrices
#' in binary format.
#'
#' @field f_path Character string specifying the file path for binary storage
#'
#' @section Methods:
#' \describe{
#'   \item{\code{new(f_path)}}{Create a new SingeCellCountData instance
#'     \itemize{
#'       \item \code{f_path}: Path to the binary file
#'     }
#'   }
#'   \item{\code{r_csc_mat_to_file(no_cells, no_genes, data, col_ptr, row_idx)}}{
#'     Write CSC matrix data to binary file
#'     \itemize{
#'       \item \code{no_cells}: Number of cells in the matrix
#'       \item \code{no_genes}: Number of genes in the matrix
#'       \item \code{data}: Vector of non-zero values
#'       \item \code{col_ptr}: Column pointer array for CSC format
#'       \item \code{row_idx}: Row indices for non-zero values
#'     }
#'   }
#'   \item{\code{file_to_r_csc_mat()}}{Read CSC matrix from binary file
#'     \itemize{
#'       \item Returns: List containing col_ptr, row_idx, data, no_cells, no_genes
#'     }
#'   }
#' }
#'
#' @name SingeCellCountData
#' @export
NULL
