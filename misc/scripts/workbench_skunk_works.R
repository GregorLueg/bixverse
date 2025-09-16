# sparse data format -----------------------------------------------------------

library(magrittr)

rextendr::document()

## csr (cell-centric) ----------------------------------------------------------

# based on the h5ad format with cells -> rows; genes -> columns
seed <- 123L
no_genes <- 100L
no_cells <- 1000L

rs_sparse_data <- rs_synthetic_sc_data_csr(
  n_genes = no_genes,
  n_cells = no_cells,
  min_genes = 25,
  max_genes = 50,
  max_exp = 50,
  seed = seed
)

length(rs_sparse_data$data)

# rm(rs_sparse_data)

# csr_matrix <- as(csr_matrix, "RsparseMatrix")

dir <- tempdir()
f_path_cells <- file.path(dir, "cells.bin")
f_path_genes <- file.path(dir, "genes.bin")

single_cell_counts <- SingeCellCountData$new(
  f_path_cells = f_path_cells,
  f_path_genes = f_path_genes
)

# list.files(dir)

# rextendr::document()

tictoc::tic()
x <- single_cell_counts$r_csr_mat_to_file(
  no_cells = no_cells,
  no_genes = no_genes,
  data = as.integer(rs_sparse_data$data),
  row_ptr = as.integer(rs_sparse_data$indptr),
  col_idx = as.integer(rs_sparse_data$indices),
  target_size = 1e5,
  min_genes = 30L
)
tictoc::toc()

tictoc::tic()
single_cell_counts$generate_gene_based_data(min_cells = 5L)
tictoc::toc()

indices <- sort(sample(1:no_cells, 10))

tictoc::tic()
return_data <- single_cell_counts$get_cells_by_indices(
  indices = order(indices),
  assay = "norm"
)
tictoc::toc()

single_cell_counts$return_full_mat(
  assay = "raw",
  cell_based = FALSE,
  verbose = TRUE
)

file.size(f_path_cells) / 1024^2
file.size(f_path_genes) / 1024^2

# csc (gene-centric) -----------------------------------------------------------

gene_indices <- sample(1:no_genes, 10)

tictoc::tic()
return_gene_data <- single_cell_counts$get_genes_by_indices(
  indices = gene_indices,
  assay = "raw"
)
tictoc::toc()

return_gene_data$row_ptr


# h5 files ---------------------------------------------------------------------

library(duckdb)

devtools::load_all()
devtools::document()

rextendr::document()

h5_path <- "~/Downloads/ERX11148735.h5ad"

bixverse_sc <- single_cell_exp(dir_data = tempdir())

bixverse_sc <- load_h5ad(bixverse_sc, h5_path = h5_path)

obs <- get_sc_obs(
  bixverse_sc,
)

bixverse_sc[1:5L, ]

bixverse_sc[[]]

bixverse_sc@dims

rust_con <- get_sc_rust_ptr(bixverse_sc)

col_index_map <- get_sc_map(bixverse_sc)[["gene_col_idx_map"]]

return_format = "cell"
cell_indices = 1:10L
assay = "raw"
.verbose = FALSE

res <- if (return_format == "cell") {
  if (is.null(cell_indices)) {
    rust_con$return_full_mat(
      assay = assay,
      cell_based = TRUE,
      verbose = .verbose
    )
  } else {
    rust_con$get_cells_by_indices(indices = cell_indices, assay = assay)
  }
} else {
  # gene
  if (is.null(gene_indices)) {
    rust_con$return_full_mat(
      assay = assay,
      cell_based = FALSE,
      verbose = .verbose
    )
  } else {
    rust_con$get_genes_by_indices(indices = gene_indices, assay = assay)
  }
}

res$indices <- col_index_map[as.character(res$indices)]

get_sc_var(bixverse_sc)

length(bixverse_sc@index_maps$cell_map)
length(bixverse_sc@index_maps$gene_col_idx_map)

bixverse_sc@dims

count_data <- res

matrix_class <- if (return_format == "cell") "dgRMatrix" else "dgCMatrix"
index_slot <- if (return_format == "cell") "j" else "i"

sparse_mat <- with(count_data, {
  args <- list(
    p = as.integer(indptr),
    x = as.numeric(data),
    Dim = as.integer(c(no_cells, no_genes))
  )
  args[[index_slot]] <- as.integer(indices)

  do.call(new, c(matrix_class, args))
})

any(is.na(col_index_map[as.character(res$indices)]))

head(obs)

bixverse_sc[[1:25]]

bixverse_sc[[c("cell_id", "lib_size", "nnz")]]

counts <- get_sc_counts(bixverse_sc, return_format = "gene")

class(counts)

counts <- bixverse_sc[,, assay = "norm", return_format = "gene"]

class(counts)

# mtx file ---------------------------------------------------------------------

rextendr::document()

dir <- tempdir()
f_path_cells <- file.path(dir, "cells.bin")
f_path_genes <- file.path(dir, "genes.bin")

single_cell_counts <- SingeCellCountData$new(
  f_path_cells = f_path_cells,
  f_path_genes = f_path_genes
)

mtx_path <- path.expand("~/Downloads/ex053/DGE.mtx")

res <- single_cell_counts$mtx_to_file(
  mtx_path = mtx_path,
  qc_params = params_sc_min_quality()
)

genes_to_keep = single_cell_counts$generate_gene_based_data(
  qc_params = params_sc_min_quality(),
  verbose = TRUE
)
tictoc::toc()

count_data <- single_cell_counts$get_cells_by_indices(1:10, assay = "norm")

count_data

max(count_data$indices)

length(genes_to_keep$gene_mask)

table(genes_to_keep$gene_mask)


single_cell_counts$get_shape()

single_cell_counts$get_genes_by_indices(indices = 1:10, assay = "norm")

rextendr::document()
devtools::load_all()

object = single_cell_exp(dir_data = tempdir())

tictoc::tic()
object = load_mtx(
  object = object,
  mtx_path = path.expand("~/Downloads/ex053/DGE.mtx"),
  obs_path = path.expand("~/Downloads/ex053/cell_metadata.csv"),
  var_path = path.expand("~/Downloads/ex053/all_genes.csv"),
  sc_qc_param = params_sc_min_quality(
    min_unique_genes = 100L,
    min_lib_size = 250L,
    min_cells = 10L
  )
)
tictoc::toc()

object[1:10L, ]
