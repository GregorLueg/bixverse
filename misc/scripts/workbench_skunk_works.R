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

rextendr::document()

object = single_cell_exp(dir_data = tempdir())

list.files(tempdir())

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

library(biomaRt)

pryr::object_size(object)

object[[c("cell_idx", "cell_id")]]

tictoc::tic()
object[[]]
tictoc::toc()

object[[1:1500]]

tictoc::tic()
counts_gene <- object[, 1:20L, assay = "norm", return_format = "gene"]
tictoc::toc()

tictoc::tic()
counts_gene <- object[, 1:20L, assay = "norm", return_format = "cell"]
tictoc::toc()

tictoc::tic()
counts_cell <- object[,, assay = "norm", return_format = "cell"]
tictoc::toc()

list.files(tempdir())

counts_gene

class(counts_gene)

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

get_sc_var(object)


get_gene_names(object)

object <- set_cell_to_keep(object, get_cell_names(object))

object <- find_hvg(object = object)

get_hvg(object)

rextendr::document()

pryr::object_size(object)

print("Using standard SVD")

pc_results <- rs_sc_pca(
  f_path_gene = get_rust_count_gene_f_path(object),
  no_pcs = 30L,
  random_svd = FALSE,
  cell_indices = get_cells_to_keep(object),
  gene_indices = get_hvg(object),
  seed = 42L,
  verbose = TRUE
)

pc_results$scores[1:5, ]

object

print("Using randomised SVD")

pc_results_randomised <- rs_sc_pca(
  f_path_gene = get_rust_count_gene_f_path(object),
  no_pcs = 30L,
  random_svd = TRUE,
  cell_indices = get_cells_to_keep(object),
  gene_indices = get_hvg(object),
  seed = 42L,
  verbose = TRUE
)

pc_results$scores[1:5, 1:5]

pc_results$loadings[1:5, 1:5]

plot(
  pc_results_randomised$scores[, 1],
  pc_results$scores[, 1],
  xlab = "Randomised SVD",
  ylab = "SVD",
  main = "PC1"
)

plot(
  pc_results_randomised$scores[, 10],
  pc_results$scores[, 10],
  xlab = "Randomised SVD",
  ylab = "SVD",
  main = "PC10"
)

plot(
  pc_results_randomised$scores[, 25],
  pc_results$scores[, 25],
  xlab = "Randomised SVD",
  ylab = "SVD",
  main = "PC25"
)

dim(pc_results_randomised$scores)
dim(pc_results_randomised$loadings)

rextendr::document()

BiocManager::install("BiocNeighbors")

diag(rs_cor2(pc_results_randomised$scores, pc_results$scores, spearman = FALSE))

dim(pc_results$scores)

tictoc::tic()
bioc_knn = BiocNeighbors::findKNN(pc_results$scores, 10, num.threads = 10L)
tictoc::toc()

bioc_knn$index[1:10, ]

bioc_knn$distance[1:10, ]

dim(bioc_knn$index)

rextendr::document()

tictoc::tic()
test_res <- rs_sc_knn(
  embd = pc_results$scores,
  no_neighbours = 10,
  seed = 101L,
  n_trees = 100L,
  search_budget = 200L,
  verbose = TRUE,
  algorithm_type = "hnsw"
)
tictoc::toc()

setDT(test_res)

head(test_res, 25)

tictoc::tic()
test_res_2 <- rs_sc_knn(
  embd = pc_results$scores,
  no_neighbours = 10,
  seed = 101L,
  n_trees = 100L,
  search_budget = 100L,
  verbose = TRUE,
  algorithm_type = "annoy"
)
tictoc::toc()

setDT(test_res_2)


head(test_res_2, 25)

length(unique(test_res$from))

length(unique(c(test_res$from, test_res$to)))

test_res[26:35]

rextendr::document()

BiocManager::install("scran")


adj_mat <- new(
  "dgRMatrix",
  p = as.integer(test_res$indptr),
  x = test_res$data,
  j = as.integer(test_res$indices),
  Dim = c(dim(pc_results$scores)[1], dim(pc_results$scores)[1])
)

g <- igraph::graph_from_data_frame(test_res_2, directed = FALSE)

communities <- igraph::cluster_louvain(g, resolution = 0.25)

table(communities$membership)

communities <- igraph::cluster_leiden(
  g,
  objective_function = "modularity",
  n_iterations = 5L
)

table(communities$membership)

setDT(test_res_2)

rextendr::document()

tictoc::tic()
leiden <- rs_leiden_clustering(
  from = as.integer(test_res_2$from),
  to = as.integer(test_res_2$to),
  weights = test_res$weight,
  max_iterations = 5L,
  res = 1,
  seed = 42L
)
tictoc::toc()

length(unique(leiden))

summary(as.numeric(table(leiden)))

library(igraph)

tictoc::tic()


tictoc::toc()

?igraph::graph_from_adjacency_matrix()

hist(table(communities$membership))

summary(as.numeric(table(communities$membership)))

length(unique(communities$membership))


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
