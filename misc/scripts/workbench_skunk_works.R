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

obs <- get_sc_obs(bixverse_sc)
var <- get_sc_var(bixverse_sc)

counts <- bixverse_sc[,, assay = "raw", return_format = "gene"]

tictoc::tic()
counts_gene <- bixverse_sc[, 1:100L, assay = "norm", return_format = "gene"]
tictoc::toc()

tictoc::tic()
counts_gene <- bixverse_sc[, 1:100L, assay = "norm", return_format = "cell"]
tictoc::toc()

get_gene_names(bixverse_sc)

library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

G_list <- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "external_gene_name", "hgnc_symbol"),
  values = var$gene_id,
  mart = mart
)
setDT(G_list)

mt_genes <- G_list[external_gene_name %like% "MT-", ensembl_gene_id]
rps_genes <- G_list[external_gene_name %like% "^RPS", ensembl_gene_id]

gs_of_interest <- list(
  "mt_perc" = mt_genes,
  "rb_perc" = rps_genes
)

bixverse_sc <- gene_set_proportions(bixverse_sc, gs_of_interest)

cells_to_keep <- bixverse_sc[[]][mt_perc <= 0.05, cell_id]

bixverse_sc <- set_cell_to_keep(bixverse_sc, cells_to_keep)

length(get_cells_to_keep(bixverse_sc))

bixverse_sc <- find_hvg(object = bixverse_sc)

rextendr::document()

pc_results_randomised <- rs_sc_pca(
  f_path_gene = get_rust_count_gene_f_path(bixverse_sc),
  no_pcs = 30L,
  random_svd = TRUE,
  cell_indices = get_cells_to_keep(bixverse_sc),
  gene_indices = get_hvg(bixverse_sc),
  seed = 42L,
  verbose = TRUE
)

dim(pc_results_randomised$scores)

# mtx file ---------------------------------------------------------------------

## demo ------------------------------------------------------------------------

rextendr::document()

# generate the object
object = single_cell_exp(dir_data = tempdir())

# load in the mtx file

# several things are happening here...
# 1.) we scan which genes to include
# 2.) we scan which cells to include
# 3.) we load in the data and save to a binarised CSR file and generate at the
# same time the normalised counts
# 4.) we take the raw counts and normalised counts from the cells and save
# them additionally in a CSC format again to disk.

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

# show size of the object
pryr::object_size(object)

# get obs table
object[[]]

# get specific columns obs table
object[[c("cell_idx", "cell_id")]]

# subset rows
object[[1:1500]]

# get vars
get_sc_var(object)

# demo count retrieval

# gene centric version
tictoc::tic()
counts_gene <- object[, 1:20L, assay = "norm", return_format = "gene"]
tictoc::toc()

class(counts_gene)

tictoc::tic()
counts_gene <- object[, 1:20L, assay = "norm", return_format = "cell"]
tictoc::toc()

class(counts_gene)

# get gene set proportions
library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

G_list <- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "external_gene_name", "hgnc_symbol"),
  values = get_gene_names(object),
  mart = mart
)
setDT(G_list)

mt_genes <- G_list[external_gene_name %like% "MT-", ensembl_gene_id]
rps_genes <- G_list[external_gene_name %like% "^RPS", ensembl_gene_id]

gs_of_interest <- list(
  "mt_perc" = mt_genes,
  "rb_perc" = rps_genes
)

object <- gene_set_proportions(object, gs_of_interest)

object[[]]

object <- set_cell_to_keep(object, get_cell_names(object))

# identify HVGs
object <- find_hvg(object = object)

get_hvg(object)

# run PCA
object <- calculate_pca_single_cell(object, no_pcs = 30L)

# demo difference PCA vs randomised SVD

# pc_results <- rs_sc_pca(
#   f_path_gene = get_rust_count_gene_f_path(object),
#   no_pcs = 30L,
#   random_svd = FALSE,
#   cell_indices = get_cells_to_keep(object),
#   gene_indices = get_hvg(object),
#   seed = 42L,
#   verbose = TRUE
# )

# pc_results_randomised <- rs_sc_pca(
#   f_path_gene = get_rust_count_gene_f_path(object),
#   no_pcs = 30L,
#   random_svd = TRUE,
#   cell_indices = get_cells_to_keep(object),
#   gene_indices = get_hvg(object),
#   seed = 42L,
#   verbose = TRUE
# )

# plot(
#   pc_results_randomised$scores[, 1],
#   pc_results$scores[, 1],
#   xlab = "Randomised SVD",
#   ylab = "SVD",
#   main = "PC1"
# )

# plot(
#   pc_results_randomised$scores[, 10],
#   pc_results$scores[, 10],
#   xlab = "Randomised SVD",
#   ylab = "SVD",
#   main = "PC10"
# )

# plot(
#   pc_results_randomised$scores[, 25],
#   pc_results$scores[, 25],
#   xlab = "Randomised SVD",
#   ylab = "SVD",
#   main = "PC25"
# )

# diag(rs_cor2(pc_results_randomised$scores, pc_results$scores, spearman = FALSE))

get_pca_factors(object)[1:5, ]

get_pca_loadings(object)[1:5, ]

devtools::load_all()

object <- find_neigbours_single_cell(
  object,
  no_embd_to_use = 50L
)

get_snn_graph(object)

get_knn_mat(object)


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
