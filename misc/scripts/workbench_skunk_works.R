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

h5_content <- rhdf5::h5ls(
  h5_path
) %>%
  data.table::setDT()

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
  sc_mtx_io_param = params_sc_mtx_io(
    path_mtx = path.expand("~/Downloads/ex053/DGE.mtx"),
    path_obs = path.expand("~/Downloads/ex053/cell_metadata.csv"),
    path_var = path.expand("~/Downloads/ex053/all_genes.csv"),
    cells_as_rows = TRUE
  ),
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
object <- find_hvg_sc(object = object, streaming = FALSE)

get_hvg(object)

# run PCA
object <- calculate_pca_sc(object, no_pcs = 30L)

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

object <- find_neighbours_sc(
  object,
  no_embd_to_use = 50L
)


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

file_res$cell_indices

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

# test tahoe 100m h5ad ---------------------------------------------------------

rextendr::document()

devtools::load_all()

tahoe_h5ad_f_path <- path.expand(
  "~/Downloads/plate1_filt_Vevo_Tahoe100M_WServicesFrom_ParseGigalab.h5ad"
)

h5_meta <- get_h5ad_dimensions(tahoe_h5ad_f_path)


object <- single_cell_exp(dir_data = tempdir())
sc_qc_param = params_sc_min_quality()
.verbose = TRUE

rust_con <- get_sc_rust_ptr(object)

tictoc::tic()
file_res <- rust_con$h5_to_file(
  cs_type = h5_meta$type,
  h5_path = path.expand(tahoe_h5ad_f_path),
  no_cells = h5_meta$dims["obs"],
  no_genes = h5_meta$dims["var"],
  qc_params = sc_qc_param,
  verbose = .verbose
)
tictoc::toc()

tictoc::tic()
rust_con$generate_gene_based_data_streaming(
  verbose = .verbose,
  batch_size = 10000L
)
tictoc::toc()

summary(file_res$cell_indices)

length(file_res$cell_indices)

length(file_res$gene_indices)

summary(file_res$gene_indices)

devtools::document()

bixverse_sc <- single_cell_exp(dir_data = tempdir())

tictoc::tic()
bixverse_sc <- stream_h5ad(object = bixverse_sc, h5_path = tahoe_h5ad_f_path)
tictoc::toc()

obs <- get_sc_obs(bixverse_sc)
var <- get_sc_var(bixverse_sc)

table(duplicated(obs$cell_idx))

summary(obs$cell_idx)

# summary(obs$rb_perc)

# library(biomaRt)
# mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

# G_list <- getBM(
#   filters = "ensembl_gene_id",
#   attributes = c("ensembl_gene_id", "external_gene_name", "hgnc_symbol"),
#   values = get_gene_names(bixverse_sc),
#   mart = mart
# )
# setDT(G_list)

# get_gene_names(bixverse_sc)

# mt_genes <- G_list[external_gene_name %like% "MT-", ensembl_gene_id]
# rps_genes <- G_list[external_gene_name %like% "^RPL", ensembl_gene_id]

gs_of_interest <- list(
  "mt_perc" = get_gene_names(bixverse_sc)[
    get_gene_names(bixverse_sc) %like% "MT-"
  ],
  "rps_perc" = get_gene_names(bixverse_sc)[
    get_gene_names(bixverse_sc) %like% "^RPS"
  ]
)

bixverse_sc <- gene_set_proportions_sc(
  bixverse_sc,
  gs_of_interest,
  streaming = TRUE
)

cells_to_keep <- bixverse_sc[[]][mt_perc <= 0.10, cell_id]

bixverse_sc <- set_cell_to_keep(bixverse_sc, cells_to_keep)

bixverse_sc <- find_hvg_sc(object = bixverse_sc, streaming = TRUE)

bixverse_sc <- calculate_pca_sc(bixverse_sc, no_pcs = 30L)

tictoc::tic()
bixverse_sc <- find_neighbours_sc(
  bixverse_sc,
  neighbours_params = params_sc_neighbours(knn_algorithm = "hnsw")
)
tictoc::toc()

rextendr::document()

pryr::object_size(bixverse_sc)

bixverse_sc@sc_map$cells_to_keep_idx

tictoc::tic()
test_count <- bixverse_sc[2857392L, , return_format = "cell", assay = "raw"]
tictoc::toc()

get_rust_count_cell_f_path(bixverse_sc)

rextendr::document()

dge_test <- rs_calculate_dge_mann_whitney(
  f_path = get_rust_count_cell_f_path(bixverse_sc),
  cell_indices_1 = 1:10000L,
  cell_indices_2 = 10001:20000L,
  min_prop = 0.05,
  TRUE
)

summary(dge_test$lfc)

hist(dge_test$lfc)

hist(dge_test$prop1)

hist(dge_test$prop2)

hist(dge_test$z_scores)

# install.packages("devtools")
devtools::install_github("immunogenomics/presto")

tictoc::tic()
test_count <- bixverse_sc[1:200000L, , return_format = "cell", assay = "norm"]
tictoc::toc()

presto_input <- Matrix::t(test_count)
presto_output <- as(presto_input, "CsparseMatrix")

y <- factor(c(rep("A", 100000), rep("B", 100000)))

tictoc::tic()
x <- presto::wilcoxauc(presto_output, y)
tictoc::toc()

tictoc::tic()
dge_test_v2 <- rs_calculate_dge_mann_whitney(
  f_path = get_rust_count_cell_f_path(bixverse_sc),
  cell_indices_1 = 1:100000L,
  cell_indices_2 = 100001:200000L,
  min_prop = 0.0,
  alternative = "greater",
  verbose = TRUE
)
tictoc::toc()

dim(x)

rextendr::document()

tictoc::tic()
dge_test <- rs_calculate_dge_mann_whitney(
  f_path = get_rust_count_cell_f_path(bixverse_sc),
  cell_indices_1 = 1:10000L,
  cell_indices_2 = 10001:20000L,
  min_prop = 0.0,
  alternative = "less",
  verbose = TRUE
)
tictoc::toc()

tictoc::tic()
dge_test_v2 <- rs_calculate_dge_mann_whitney(
  f_path = get_rust_count_cell_f_path(bixverse_sc),
  cell_indices_1 = 1:10000L,
  cell_indices_2 = 10001:20000L,
  min_prop = 0.0,
  alternative = "greater",
  verbose = TRUE
)
tictoc::toc()

tictoc::tic()
dge_test_v3 <- rs_calculate_dge_mann_whitney(
  f_path = get_rust_count_cell_f_path(bixverse_sc),
  cell_indices_1 = 1:10000L,
  cell_indices_2 = 10001:20000L,
  min_prop = 0.0,
  alternative = "twosided",
  verbose = TRUE
)
tictoc::toc()

plot(dge_test$p_values, dge_test$p_values)

plot(dge_test_v3$p_values, dge_test_v2$p_values)

class(x)

plot(x[x$group == "A", ]$pval, dge_test_v3$p_values)

plot(x[x$group == "A", ]$logFC, dge_test$lfc)

cor(x[x$group == "A", ]$logFC, dge_test$lfc, method = "spearman")

plot(x[x$group == "A", ]$logFC[1:10], dge_test$lfc[1:10])

x[x$group == "A", ]$pval[1:10]

dge_test$p_values[1:10]

gc()

?wilcox.test

# debug function ---------------------------------------------------------------

devtools::load_all()

get_h5ad_dimensions(tahoe_h5ad_f_path)

get_h5ad_dimensions(h5_path)

rextendr::document()

h5_content <- rhdf5::h5ls(
  h5_path
) %>%
  data.table::setDT()

no_obs <- h5_content[
  group == "/obs" & otype == "H5I_DATASET"
][1, as.numeric(dim)]


no_var <- h5_content[
  group == "/var" & otype == "H5I_DATASET"
][1, as.numeric(dim)]

indptr <- h5_content[
  group == "/X" & name == "indptr",
  as.numeric(dim)
]

cs_format <- ifelse(no_var + 1 == indptr, "CSR", "CSC")

# version 2

h5_content <- rhdf5::h5ls(
  tahoe_h5ad_f_path
) %>%
  data.table::setDT()

h5_path <- "~/Downloads/ERX11148735.h5ad"

h5_content_2 <- rhdf5::h5ls(
  h5_path
) %>%
  data.table::setDT()


no_var <- h5_content[
  group == "/var" & otype == "H5I_DATASET"
][1, as.numeric(dim)]

indptr <- h5_content[
  group == "/X" & name == "indptr",
  as.numeric(dim)
]


# generating test data ---------------------------------------------------------

rextendr::document()

n_cells = 1000L
n_genes = 100L
n_background_genes_exp = c(3L, 15L)
background_exp_range = c(5L, 10L)

marker_genes <- list(
  cell_type_1 = list(
    marker_genes = 0:9L,
    marker_exp_range = c(10L, 50L),
    markers_per_cell = c(2L, 8L)
  ),
  cell_type_2 = list(
    marker_genes = 10:19L,
    marker_exp_range = c(10L, 50L),
    markers_per_cell = c(2L, 8L)
  ),
  cell_type_3 = list(
    marker_genes = 20:29L,
    marker_exp_range = c(10L, 50L),
    markers_per_cell = c(2L, 8L)
  )
)

data <- rs_synthetic_sc_data_with_cell_types(
  n_cells = n_cells,
  n_genes = n_genes,
  n_background_genes_exp = n_background_genes_exp,
  background_exp_range = n_background_genes_exp,
  cell_configs = marker_genes,
  seed = 42L
)

counts <- new(
  "dgRMatrix",
  p = as.integer(data$indptr),
  x = as.numeric(data$data),
  j = as.integer(data$indices),
  Dim = as.integer(c(n_cells, n_genes))
)

rownames(counts) <- sprintf("cell_%04d", 1:1000)
colnames(counts) <- sprintf("gene_%03d", 1:100)

obs <- data.table(
  cell_id = sprintf("cell_%04d", 1:1000),
  cell_grp = sprintf("cell_type_%i", data$cell_type_indices + 1)
)

var <- data.table(
  gene_id = sprintf("gene_%03d", 1:100),
  ensembl_id = sprintf("ens_%03d", 1:100)
)


obs[[1]]

class(counts)

single_cell_test_data <- generate_single_cell_test_data()

devtools::load_all()

f_path = file.path(tempdir(), "csr_test.h5ad")

write_h5ad_sc(
  f_path = f_path,
  counts = single_cell_test_data$counts,
  obs = single_cell_test_data$obs,
  single_cell_test_data$var,
  .verbose = FALSE
)
