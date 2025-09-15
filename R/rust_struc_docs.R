# #' @name SingeCellCountData
# #' @title Single Cell Count Data Handler
# #'
# #' @description
# #' A class for efficiently reading and writing single cell count matrices
# #' in binary format. Supports both raw counts and normalized data with
# #' optimized disk I/O for large-scale single cell datasets.
# #'
# #' @section Methods:
# #' \subsection{Public methods}{
# #' \itemize{
# #' \item \href{#method-SingeCellCountData-new}{\code{SingeCellCountData$new()}}
# #' \item \href{#method-SingeCellCountData-r_csc_mat_to_file}{\code{SingeCellCountData$r_csc_mat_to_file()}}
# #' \item \href{#method-SingeCellCountData-file_to_r_csc_mat}{\code{SingeCellCountData$file_to_r_csc_mat()}}
# #' \item \href{#method-SingeCellCountData-get_cells_by_indices}{\code{SingeCellCountData$get_cells_by_indices()}}
# #' }
# #' }
# #' \if{html}{\out{<hr>}}
# #' \if{html}{\out{<a id="method-SingeCellCountData-new"></a>}}
# #' \if{latex}{\out{\hypertarget{method-SingeCellCountData-new}{}}}
# #' \subsection{Method \code{new()}}{
# #' Create a new SingeCellCountData instance
# #' \subsection{Usage}{
# #' \if{html}{\out{<div class="r">}}\preformatted{SingeCellCountData$new(f_path_cells, f_path_genes)}\if{html}{\out{</div>}}
# #' }
# #'
# #' \subsection{Arguments}{
# #' \if{html}{\out{<div class="arguments">}}
# #' \describe{
# #' \item{\code{f_path_cells}}{Character string. Path to the binary file for cell data storage.}
# #' \item{\code{f_path_genes}}{Character string. Path to the binary file for gene data storage.}
# #' }
# #' \if{html}{\out{</div>}}
# #' }
# #' \subsection{Returns}{
# #' Returns the initialized SingeCellCountData instance.
# #' }
# #' }
# #' \if{html}{\out{<hr>}}
# #' \if{html}{\out{<a id="method-SingeCellCountData-r_csc_mat_to_file"></a>}}
# #' \if{latex}{\out{\hypertarget{method-SingeCellCountData-r_csc_mat_to_file}{}}}
# #' \subsection{Method \code{r_csc_mat_to_file()}}{
# #' Write CSC matrix data to binary file with normalization
# #' \subsection{Usage}{
# #' \if{html}{\out{<div class="r">}}\preformatted{SingeCellCountData$r_csc_mat_to_file(no_cells, no_genes, data, col_ptr, row_idx, target_size)}\if{html}{\out{</div>}}
# #' }
# #'
# #' \subsection{Arguments}{
# #' \if{html}{\out{<div class="arguments">}}
# #' \describe{
# #' \item{\code{no_cells}}{Integer. Number of cells (columns) in the matrix.}
# #' \item{\code{no_genes}}{Integer. Number of genes (rows) in the matrix.}
# #' \item{\code{data}}{Integer vector. Non-zero values from the CSC matrix.}
# #' \item{\code{col_ptr}}{Integer vector. Column pointer array for CSC format.}
# #' \item{\code{row_idx}}{Integer vector. Row indices for non-zero values.}
# #' \item{\code{target_size}}{Numeric. Target library size for normalization (e.g., 1e6 for CPM).}
# #' }
# #' \if{html}{\out{</div>}}
# #' }
# #' }
# #' \if{html}{\out{<hr>}}
# #' \if{html}{\out{<a id="method-SingeCellCountData-file_to_r_csc_mat"></a>}}
# #' \if{latex}{\out{\hypertarget{method-SingeCellCountData-file_to_r_csc_mat}{}}}
# #' \subsection{Method \code{file_to_r_csc_mat()}}{
# #' Read complete CSC matrix from binary file
# #' \subsection{Usage}{
# #' \if{html}{\out{<div class="r">}}\preformatted{SingeCellCountData$file_to_r_csc_mat(assay)}\if{html}{\out{</div>}}
# #' }
# #'
# #' \subsection{Arguments}{
# #' \if{html}{\out{<div class="arguments">}}
# #' \describe{
# #' \item{\code{assay}}{Character string. Data type to return: \code{"raw"} for raw counts or \code{"norm"} for normalized data.}
# #' }
# #' \if{html}{\out{</div>}}
# #' }
# #' \subsection{Returns}{
# #' List containing:
# #' \itemize{
# #' \item \code{col_ptr}: Column pointer vector for CSC format
# #' \item \code{row_idx}: Row indices of non-zero values
# #' \item \code{data}: Values vector (integer for raw, numeric for normalized)
# #' \item \code{no_cells}: Number of cells
# #' \item \code{no_genes}: Number of genes
# #' }
# #' }
# #' }
# #' \if{html}{\out{<hr>}}
# #' \if{html}{\out{<a id="method-SingeCellCountData-get_cells_by_indices"></a>}}
# #' \if{latex}{\out{\hypertarget{method-SingeCellCountData-get_cells_by_indices}{}}}
# #' \subsection{Method \code{get_cells_by_indices()}}{
# #' Read specific cells by their indices
# #' \subsection{Usage}{
# #' \if{html}{\out{<div class="r">}}\preformatted{SingeCellCountData$get_cells_by_indices(indices, assay)}\if{html}{\out{</div>}}
# #' }
# #'
# #' \subsection{Arguments}{
# #' \if{html}{\out{<div class="arguments">}}
# #' \describe{
# #' \item{\code{indices}}{Integer vector. Cell indices to retrieve (1-based indexing).}
# #' \item{\code{assay}}{Character string. Data type to return: \code{"raw"} for raw counts or \code{"norm"} for normalized data.}
# #' }
# #' \if{html}{\out{</div>}}
# #' }
# #' \subsection{Returns}{
# #' List containing:
# #' \itemize{
# #' \item \code{col_ptr}: Column pointer vector for selected cells
# #' \item \code{row_idx}: Row indices of non-zero values
# #' \item \code{data}: Values vector (integer for raw, numeric for normalized)
# #' \item \code{no_cells}: Number of selected cells
# #' \item \code{no_genes}: Total number of genes
# #' }
# #' }
# #' }
# #'
# #' @export
# NULL
