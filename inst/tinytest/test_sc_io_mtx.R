# mtx io -----------------------------------------------------------------------

test_temp_dir <- file.path(
  tempdir(),
  paste0("test_", format(Sys.time(), "%Y%m%d_%H%M%S_"), sample(1000:9999, 1))
)
dir.create(test_temp_dir, recursive = TRUE)

## parameters ------------------------------------------------------------------

# testing parameters
min_lib_size <- 300L
min_genes_exp <- 45L
min_cells_exp <- 500L

## synthetic data --------------------------------------------------------------

single_cell_test_data <- generate_single_cell_test_data()

f_path_v1 <- file.path(test_temp_dir, "cells_csv")
f_path_v2 <- file.path(test_temp_dir, "genes_tsv")

dir.create(f_path_v1, showWarnings = FALSE, recursive = TRUE)
dir.create(f_path_v2, showWarnings = FALSE, recursive = TRUE)

counts_csc <- as(single_cell_test_data$counts, "CsparseMatrix")

# save a version with rows = cells and format type csv for the rest

write_cellranger_output(
  f_path = f_path_v1,
  counts = single_cell_test_data$counts,
  obs = single_cell_test_data$obs,
  var = single_cell_test_data$var,
  rows = "cells",
  format_type = "csv",
  .verbose = FALSE
)

# save a version with rows = genes and format type tsv for the rest

write_cellranger_output(
  f_path = f_path_v2,
  counts = single_cell_test_data$counts,
  obs = single_cell_test_data$obs,
  var = single_cell_test_data$var,
  rows = "genes",
  format_type = "tsv",
  .verbose = FALSE
)

### expected data --------------------------------------------------------------

genes_pass <- which(
  Matrix::colSums(single_cell_test_data$counts != 0) >= min_cells_exp
)

cells_pass <- which(
  (Matrix::rowSums(single_cell_test_data$counts[, genes_pass]) >=
    min_lib_size) &
    (Matrix::rowSums(single_cell_test_data$counts[, genes_pass] != 0) >=
      min_genes_exp)
)

expect_true(
  current = length(genes_pass) > 80 & length(genes_pass) != 100,
  info = "mtx - sensible amount of genes pass"
)

expect_true(
  current = length(cells_pass) > 800 & length(cells_pass) != 1000,
  info = "mtx - sensible amount of cells pass"
)

counts_filtered <- single_cell_test_data$counts[cells_pass, genes_pass]

counts_filtered_csc <- counts_csc[cells_pass, genes_pass]

obs_filtered <- single_cell_test_data$obs[cells_pass, ]

vars_filtered <- single_cell_test_data$var[genes_pass, ]

sc_qc_param = params_sc_min_quality(
  min_unique_genes = min_genes_exp,
  min_lib_size = min_lib_size,
  min_cells = min_cells_exp
)

params_cells_rows_csv <- params_sc_mtx_io(
  path_mtx = file.path(f_path_v1, "mat.mtx"),
  path_obs = file.path(f_path_v1, "barcodes.csv"),
  path_var = file.path(f_path_v1, "features.csv"),
  cells_as_rows = TRUE,
  has_hdr = TRUE
)

params_genes_rows_tsv <- params_sc_mtx_io(
  path_mtx = file.path(f_path_v2, "mat.mtx"),
  path_obs = file.path(f_path_v2, "barcodes.tsv"),
  path_var = file.path(f_path_v2, "features.tsv"),
  cells_as_rows = FALSE,
  has_hdr = TRUE
)

### param checks ---------------------------------------------------------------

expect_true(
  current = bixverse:::checkScMtxIO(params_cells_rows_csv),
  info = "MTX wrapper tests work as anticipated"
)

expect_true(
  current = bixverse:::checkScMtxIO(params_genes_rows_tsv),
  info = "MTX wrapper tests work as anticipated (v2)"
)

### row = cells ; csv format ---------------------------------------------------

#### rust ----------------------------------------------------------------------

# test the underlying rust directly
sc_object <- suppressWarnings(single_cell_exp(dir_data = test_temp_dir))

rust_con <- get_sc_rust_ptr(sc_object)

file_res <- with(
  params_cells_rows_csv,
  rust_con$mtx_to_file(
    mtx_path = path_mtx,
    qc_params = sc_qc_param,
    cells_as_rows = cells_as_rows,
    verbose = FALSE
  )
)

expect_equivalent(
  current = file_res$cell_indices + 1,
  target = cells_pass,
  info = paste("mtx to binary - correct cells being kept")
)

expect_equivalent(
  current = file_res$gene_indices + 1,
  target = genes_pass,
  info = paste("mtx to binary - correct genes being kept")
)

#### duckdb csv ----------------------------------------------------------------

duckdb_con <- get_sc_duckdb(sc_object)

with(
  params_cells_rows_csv,
  {
    duckdb_con$populate_obs_from_plain_text(
      f_path = path_obs,
      filter = as.integer(file_res$cell_indices + 1),
      has_hdr = TRUE
    )
    duckdb_con$populate_var_from_plain_text(
      f_path = path_var,
      filter = as.integer(file_res$gene_indices + 1),
      has_hdr = TRUE
    )
  }
)

obs_filtered_db <- duckdb_con$get_obs_table()

vars_filtered_db <- duckdb_con$get_vars_table()

expect_equivalent(
  current = obs_filtered_db$cell_id,
  target = obs_filtered$cell_id,
  info = "obs table from CSV - correct cells kept"
)

expect_equivalent(
  current = vars_filtered_db$gene_id,
  target = vars_filtered$gene_id,
  info = "obs table from CSV - correct genes kept"
)

#### full object ---------------------------------------------------------------

sc_object <- suppressWarnings(single_cell_exp(dir_data = test_temp_dir))

sc_object <- load_mtx(
  object = sc_object,
  sc_mtx_io_param = params_cells_rows_csv,
  sc_qc_param = sc_qc_param,
  streaming = FALSE,
  .verbose = FALSE
)

##### observations -------------------------------------------------------------

obs_table_obj <- sc_object[[]]

expect_equal(
  current = obs_table_obj$cell_id,
  target = obs_filtered$cell_id,
  info = "obs table from csv (from object) - correct cells kept"
)

expect_equal(
  current = obs_table_obj$cell_grp,
  target = obs_filtered$cell_grp,
  info = "obs table from csv (from object) - cell group correct"
)

expect_equivalent(
  current = obs_table_obj$lib_size,
  target = Matrix::rowSums(counts_filtered),
  info = "obs table from csv (from object) - library size correct"
)

expect_equivalent(
  current = obs_table_obj$nnz,
  target = Matrix::rowSums(counts_filtered != 0),
  info = "obs table from csv (from object) - nnz correct"
)

#### vars ----------------------------------------------------------------------

vars_object <- get_sc_var(sc_object)

expect_equal(
  current = vars_object$gene_id,
  target = vars_filtered$gene_id,
  info = "var table from CSV (from object) - correct genes kept"
)

#### counts --------------------------------------------------------------------

counts_csr_obj <- sc_object[,, return_format = "cell"]

expect_equal(
  current = counts_csr_obj,
  target = counts_filtered,
  info = paste("count retrieval - full CSR")
)

set.seed(42L)

random_cell_indices <- sample(1:nrow(counts_csr_obj), size = 20L)

expect_equal(
  current = sc_object[random_cell_indices, , return_format = "cell"],
  target = counts_filtered[random_cell_indices, ],
  info = paste("count retrieval - CSR - random index selection")
)

counts_csc_object <- sc_object[,, return_format = "gene"]

expect_equal(
  current = counts_csc_object,
  target = counts_filtered_csc,
  info = paste("count retrieval - full CSC")
)

set.seed(123L)

random_gene_indices <- sample(1:ncol(counts_csr_obj), size = 20L)

expect_equal(
  current = sc_object[, random_gene_indices, return_format = "gene"],
  target = counts_filtered_csc[, random_gene_indices],
  info = paste("count retrieval - full CSC - random index selection")
)

#### test getters --------------------------------------------------------------

cell_names_obj <- get_cell_names(sc_object)
gene_names_obj <- get_gene_names(sc_object)

expect_equal(
  current = cell_names_obj,
  target = obs_filtered$cell_id,
  info = "correct cell names"
)

expect_equal(
  current = gene_names_obj,
  target = vars_filtered$gene_id,
  info = "correct gene names"
)

### rust = genes ; tsv format --------------------------------------------------

sc_object <- suppressWarnings(single_cell_exp(dir_data = test_temp_dir))

#### rust ----------------------------------------------------------------------

rust_con <- get_sc_rust_ptr(sc_object)

file_res <- with(
  params_genes_rows_tsv,
  rust_con$mtx_to_file(
    mtx_path = path_mtx,
    qc_params = sc_qc_param,
    cells_as_rows = cells_as_rows,
    verbose = FALSE
  )
)

expect_equivalent(
  current = file_res$cell_indices + 1,
  target = cells_pass,
  info = paste("mtx to binary - correct cells being kept")
)

expect_equivalent(
  current = file_res$gene_indices + 1,
  target = genes_pass,
  info = paste("mtx to binary - correct genes being kept")
)

#### duckdb tsv ----------------------------------------------------------------

duckdb_con <- get_sc_duckdb(sc_object)

with(
  params_genes_rows_tsv,
  {
    duckdb_con$populate_obs_from_plain_text(
      f_path = path_obs,
      filter = as.integer(file_res$cell_indices + 1),
      has_hdr = TRUE
    )
    duckdb_con$populate_var_from_plain_text(
      f_path = path_var,
      filter = as.integer(file_res$gene_indices + 1),
      has_hdr = TRUE
    )
  }
)

obs_filtered_db <- duckdb_con$get_obs_table()

vars_filtered_db <- duckdb_con$get_vars_table()

expect_equivalent(
  current = obs_filtered_db$cell_id,
  target = obs_filtered$cell_id,
  info = "obs table from TSV - correct cells kept"
)

expect_equivalent(
  current = vars_filtered_db$gene_id,
  target = vars_filtered$gene_id,
  info = "obs table from TSV - correct genes kept"
)

#### full object ---------------------------------------------------------------

sc_object <- suppressWarnings(single_cell_exp(dir_data = test_temp_dir))

sc_object <- load_mtx(
  object = sc_object,
  sc_mtx_io_param = params_genes_rows_tsv,
  sc_qc_param = sc_qc_param,
  streaming = FALSE,
  .verbose = FALSE
)

##### observations -------------------------------------------------------------

obs_table_obj <- sc_object[[]]

expect_equal(
  current = obs_table_obj$cell_id,
  target = obs_filtered$cell_id,
  info = "obs table from csv (from object) - correct cells kept"
)

expect_equal(
  current = obs_table_obj$cell_grp,
  target = obs_filtered$cell_grp,
  info = "obs table from csv (from object) - cell group correct"
)

expect_equivalent(
  current = obs_table_obj$lib_size,
  target = Matrix::rowSums(counts_filtered),
  info = "obs table from csv (from object) - library size correct"
)

expect_equivalent(
  current = obs_table_obj$nnz,
  target = Matrix::rowSums(counts_filtered != 0),
  info = "obs table from csv (from object) - nnz correct"
)

#### vars ----------------------------------------------------------------------

vars_object <- get_sc_var(sc_object)

expect_equal(
  current = vars_object$gene_id,
  target = vars_filtered$gene_id,
  info = "var table from CSV (from object) - correct genes kept"
)

#### counts --------------------------------------------------------------------

counts_csr_obj <- sc_object[,, return_format = "cell"]

expect_equal(
  current = counts_csr_obj,
  target = counts_filtered,
  info = paste("count retrieval - full CSR")
)

set.seed(42L)

random_cell_indices <- sample(1:nrow(counts_csr_obj), size = 20L)

expect_equal(
  current = sc_object[random_cell_indices, , return_format = "cell"],
  target = counts_filtered[random_cell_indices, ],
  info = paste("count retrieval - CSR - random index selection")
)

counts_csc_object <- sc_object[,, return_format = "gene"]

expect_equal(
  current = counts_csc_object,
  target = counts_filtered_csc,
  info = paste("count retrieval - full CSC")
)

set.seed(123L)

random_gene_indices <- sample(1:ncol(counts_csr_obj), size = 20L)

expect_equal(
  current = sc_object[, random_gene_indices, return_format = "gene"],
  target = counts_filtered_csc[, random_gene_indices],
  info = paste("count retrieval - full CSC - random index selection")
)

#### test getters --------------------------------------------------------------

cell_names_obj <- get_cell_names(sc_object)
gene_names_obj <- get_gene_names(sc_object)

expect_equal(
  current = cell_names_obj,
  target = obs_filtered$cell_id,
  info = "correct cell names"
)

expect_equal(
  current = gene_names_obj,
  target = vars_filtered$gene_id,
  info = "correct gene names"
)

on.exit(unlink(test_temp_dir, recursive = TRUE, force = TRUE), add = TRUE)
