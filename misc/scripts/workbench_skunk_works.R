# sparse data format -----------------------------------------------------------

library(magrittr)

rextendr::document()

## csc (cell-centric) ----------------------------------------------------------

seed <- 123L
no_genes <- 20000L
no_cells <- 500000L

rs_sparse_data <- rs_synthetic_sc_data(
  n_genes = no_genes,
  n_cells = no_cells,
  min_genes = 200,
  max_genes = 1000,
  max_exp = 50,
  seed = seed
)

csc_matrix <- Matrix::sparseMatrix(
  i = rs_sparse_data$row_indices + 1,
  p = rs_sparse_data$col_ptrs,
  x = rs_sparse_data$data,
  dims = c(no_genes, no_cells),
  dimnames = list(
    sprintf("gene_%i", 1:no_genes),
    sprintf("cell_%i", 1:no_cells)
  )
)

dim(csc_matrix)

# check sparsity
(length(csc_matrix@x) / 1000) / ((no_genes / 1000) * (no_cells / 1000))

dir <- tempdir()
f_path_cells <- file.path(dir, "cells.bin")
f_path_genes <- file.path(dir, "genes.bin")

single_cell_counts <- SingeCellCountData$new(
  f_path_cells = f_path_cells,
  f_path_genes = f_path_genes
)

tictoc::tic()
single_cell_counts$r_csc_mat_to_file(
  no_cells = csc_matrix@Dim[2],
  no_genes = csc_matrix@Dim[1],
  data = as.integer(csc_matrix@x),
  col_ptr = csc_matrix@p,
  row_idx = csc_matrix@i,
  target_size = 1e5
)
tictoc::toc()

file.info(f_path_cells)$size / 1024^2

tictoc::tic()
return_data <- single_cell_counts$file_to_r_csc_mat(assay = "raw")
tictoc::toc()

indices <- sample(1:no_cells, 10000)

tictoc::tic()
return_data <- single_cell_counts$get_cells_by_indices(
  indices = indices,
  assay = "norm"
)
tictoc::toc()

# csr (gene-centric) -----------------------------------------------------------

## TODO

csr_matrix <- as(
  Matrix::Matrix(count_mat, sparse = TRUE),
  "RsparseMatrix"
)

length(csc_matrix@p)

csr_matrix@p
csr_matrix@j
csr_matrix@x


rextendr::document()
