# # sparse data structures -------------------------------------------------------

# ## csr data --------------------------------------------------------------------

# no_cells <- 1000
# no_genes <- 500

# rs_sparse_data <- rs_synthetic_sc_data_csr(
#   n_genes = no_genes,
#   n_cells = no_cells,
#   min_genes = 20,
#   max_genes = 50,
#   max_exp = 50,
#   seed = 123L
# )

# csr_matrix <- Matrix::sparseMatrix(
#   j = rs_sparse_data$indices + 1,
#   p = rs_sparse_data$indptr,
#   x = rs_sparse_data$data,
#   dims = c(no_cells, no_genes),
#   dimnames = list(
#     sprintf("cell_%i", 1:no_cells),
#     sprintf("gene_%i", 1:no_genes)
#   )
# )
# csr_matrix <- as(csr_matrix, "RsparseMatrix")

# ## write file ------------------------------------------------------------------

# dir <- tempdir()
# f_path_cells <- file.path(dir, "cell_data.bin")
# f_path_genes <- file.path(dir, "gene_data.bin")

# single_cell_counts <- SingeCellCountData$new(
#   f_path_cells = f_path_cells,
#   f_path_genes = f_path_genes
# )

# single_cell_counts$r_csr_mat_to_file(
#   no_cells = csr_matrix@Dim[1],
#   no_genes = csr_matrix@Dim[2],
#   data = as.integer(csr_matrix@x),
#   row_ptr = csr_matrix@p,
#   col_idx = csr_matrix@j,
#   target_size = 1e5
# )

# expect_true(
#   current = "cell_data.bin" %in% list.files(dir),
#   info = "test file for cells was written"
# )

# ## get full data back ----------------------------------------------------------

# full_data <- single_cell_counts$file_to_r_csr_mat(assay = "raw")

# expect_equal(
#   current = full_data$row_ptr,
#   target = csr_matrix@p,
#   info = "CSR format - Recovered row pointers"
# )

# expect_equal(
#   current = full_data$col_idx,
#   target = csr_matrix@j,
#   info = "CSR format - Recovered col indices"
# )

# expect_equal(
#   current = full_data$data,
#   target = csr_matrix@x,
#   info = "CSR format - Recovered data"
# )

# ## get cells by indices --------------------------------------------------------

# indices <- c(1, 2, 10, 15)

# r_subsetted_matrix <- csr_matrix[indices, ]

# rust_subsetted_data <- single_cell_counts$get_cells_by_indices(
#   indices = as.integer(indices),
#   assay = "raw"
# )

# expect_equal(
#   current = rust_subsetted_data$row_ptr,
#   target = r_subsetted_matrix@p,
#   info = "CSR format - Recovered column pointers - index subsetting"
# )

# expect_equal(
#   current = rust_subsetted_data$col_idx,
#   target = r_subsetted_matrix@j,
#   info = "CSR format - Recovered row indices - subset"
# )

# expect_equal(
#   current = rust_subsetted_data$data,
#   target = r_subsetted_matrix@x,
#   info = "CSR format - Recovered data - subset"
# )

# ## generate the gene file ------------------------------------------------------

# single_cell_counts$generate_gene_based_data()

# expect_true(
#   current = "gene_data.bin" %in% list.files(dir),
#   info = "test file for cells was written"
# )

# ## test subsetting in the gene file --------------------------------------------

# gene_indices <- c(2, 4, 6, 8, 100)

# # r version
# csc_matrix <- csr_matrix[, gene_indices]
# csc_matrix <- as(csc_matrix, "CsparseMatrix")

# rust_subsetted_data_genes <- single_cell_counts$get_genes_by_indices(
#   indices = as.integer(gene_indices),
#   assay = "raw"
# )

# expect_equal(
#   current = rust_subsetted_data_genes$row_ptr,
#   target = csc_matrix@p,
#   info = "CSC format - Recovered row pointers - gene index subsetting"
# )

# expect_equal(
#   current = rust_subsetted_data_genes$col_idx,
#   target = csc_matrix@i,
#   info = "CSC format - Recovered column indices - gene index subsetting"
# )

# expect_equal(
#   current = rust_subsetted_data_genes$data,
#   target = csc_matrix@x,
#   info = "CSC format - Recovered data - gene index subsetting"
# )
