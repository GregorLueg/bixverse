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

# cells x genes

bixverse_sc <- single_cell_exp(dir_data = tempdir())

bixverse_sc <- load_h5ad(object = bixverse_sc, h5_path = h5_path)


S7::method(load_h5ad, single_cell_exp)(new_name, h5_path = h5_path)

obs <- get_sc_obs(
  bixverse_sc
)

get_sc_map(bixverse_sc)

var <- get_sc_var(bixverse_sc)

counts <- bixverse_sc[,, assay = "raw", return_format = "gene"]

tictoc::tic()
counts_gene <- bixverse_sc[, 1:100L, assay = "norm", return_format = "gene"]
tictoc::toc()

tictoc::tic()
counts_gene <- bixverse_sc[, 1:100L, assay = "norm", return_format = "cell"]
tictoc::toc()

var <- get_sc_var(bixverse_sc)

library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

G_list <- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "external_gene_name"),
  values = var$gene_id,
  mart = mart
)

setDT(G_list)

mt_genes <- G_list[external_gene_name %like% "MT-", ensembl_gene_id]

lib_size <- Matrix::rowSums(counts)

mt_reads <- Matrix::rowSums(counts[, mt_genes])

r_mt <- mt_reads / lib_size

x <- bixverse_sc@sc_map

r_gene_idx <- x[["gene_mapping"]][gene_ids]

x$gene_column_index_map[as.character(r_gene_idx)]

rust_gene_idx <- as.integer(names(x$gene_column_index_map[as.character(
  r_gene_idx
)]))

get_updated_gene_indices(x, rust_gene_idx)

r_gene_idx

bixverse_sc@sc_map$gene_column_index_map

bixverse_sc@sc_map$gene_mapping

gene_idx_of_interest <- as.numeric(names(bixverse_sc@sc_map$gene_column_index_map[bixverse_sc@sc_map$gene_mapping[
  mt_genes
]]))

gs_of_interest <- list("mt_perc" = as.integer(gene_idx_of_interest))

rextendr::document()

tictoc::tic()
qc_res <- rs_sc_get_gene_set_perc(
  f_path_cell = file.path(bixverse_sc@dir_data, "counts_cells.bin"),
  gene_set_idx = gs_of_interest
)
tictoc::toc()

qc_res$mt_perc$percentage[1:10]
r_mt[1:10]

qc_res$mt_perc$sum[1:10]
mt_reads[1:10]

qc_res$mt_perc$lib_size[1:10]
lib_size[1:10]

devtools::load_all()

bixverse_sc[["mt_perc"]] <- qc_res$mt_perc

bixverse_sc[[]]

bixverse_sc[[]][mt_perc <= 0.05]

tictoc::tic()
obs <- get_sc_obs(
  bixverse_sc
)
tictoc::toc()


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

# data <- (idxptr, idx, val, val)

# Seurat (EVERYTHING is double precision - f64)
# raw counts -> (idxptr, idx, val)
# norm counts -> (idxptr, idx, val)
# SCT counts -> (idxptr, idx, val)
# SCT data -> (idxptr, idx, val)
# scaled data (dense) -> 2500 HVGs -> N x n_cells (dense!)

# Rebuild of the h5ad parsing --------------------------------------------------

rextendr::document()
devtools::load_all()

h5_path <- path.expand("~/Downloads/ERX11148735.h5ad")
h5_meta <- get_h5ad_dimensions(h5_path)

object <- single_cell_exp(dir_data = tempdir())
h5_path
sc_qc_param = params_sc_min_quality()
.verbose = TRUE

h5_meta <- get_h5ad_dimensions(h5_path)

rust_con <- get_sc_rust_ptr(object)

file_res <- rust_con$h5_to_file(
  cs_type = h5_meta$type,
  h5_path = path.expand(h5_path),
  no_cells = h5_meta$dims["obs"],
  no_genes = h5_meta$dims["var"],
  qc_params = sc_qc_param,
  verbose = .verbose
)

rust_con$generate_gene_based_data(
  qc_params = sc_qc_param,
  verbose = .verbose
)

duckdb_con <- get_sc_duckdb(object)
duckdb_con$populate_obs_from_h5(
  h5_path = h5_path,
  filter = as.integer(file_res$cell_indices + 1)
)
duckdb_con$populate_vars_from_h5(
  h5_path = h5_path,
  filter = as.integer(file_res$gene_indices + 1)
)

length(file_res$cell_indices)

duckdb_con$get_obs_table()

cell_res_dt <- data.table::setDT(file_res[c("nnz", "lib_size")])

duckdb_con$add_data_obs(new_data = cell_res_dt)
cell_map <- duckdb_con$get_obs_index_map()
gene_map <- duckdb_con$get_var_index_map()

S7::prop(object, "dims") <- as.integer(rust_con$get_shape())
object <- set_cell_mapping(x = object, cell_map = cell_map)
object <- set_gene_mapping(x = object, gene_map = gene_map)


# Bebuild the mtx parsing ------------------------------------------------------

rextendr::document()

dir <- tempdir()
f_path_cells <- file.path(dir, "cells.bin")
f_path_genes <- file.path(dir, "genes.bin")

single_cell_counts <- SingeCellCountData$new(
  f_path_cells = f_path_cells,
  f_path_genes = f_path_genes
)

mtx_path <- path.expand("~/Downloads/ex053/DGE.mtx")

tictoc::tic()
res <- single_cell_counts$mtx_to_file(
  mtx_path = mtx_path,
  qc_params = params_sc_min_quality(),
  verbose = TRUE
)
tictoc::toc()

gene_res <- single_cell_counts$generate_gene_based_data(
  verbose = TRUE
)

single_cell_counts$get_cells_by_indices(indices = 1:10L, assay = "norm")

single_cell_counts$get_genes_by_indices(indices = 1:25L, assay = "raw")

length(res$cell_indices)

length(res$lib_size)
