# h5 io ------------------------------------------------------------------------

library(magrittr)

test_temp_dir <- file.path(
  tempdir(),
  "io_h5"
)

dir.create(test_temp_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot("Test directory does not exist" = dir.exists(test_temp_dir))

## parameters ------------------------------------------------------------------

# testing parameters
min_lib_size <- 300L
min_genes_exp <- 45L
min_cells_exp <- 500L

## synthetic data --------------------------------------------------------------

single_cell_test_data <- generate_single_cell_test_data()

# CSR version
f_path_csr = file.path(test_temp_dir, "csr_test.h5ad")

write_h5ad_sc(
  f_path = f_path_csr,
  counts = single_cell_test_data$counts,
  obs = single_cell_test_data$obs,
  single_cell_test_data$var,
  .verbose = FALSE
)

# CSC version
counts_csc <- as(single_cell_test_data$counts, "CsparseMatrix")

f_path_csc = file.path(test_temp_dir, "csc_test.h5ad")

write_h5ad_sc(
  f_path = f_path_csc,
  counts = counts_csc,
  obs = single_cell_test_data$obs,
  single_cell_test_data$var,
  .verbose = FALSE
)

# tests ------------------------------------------------------------------------

## format identification -------------------------------------------------------

expected_dims <- setNames(c(1000, 100), c("obs", "var"))

h5_meta_csr <- bixverse:::get_h5ad_dimensions(f_path_csr)

h5_meta_csc <- bixverse:::get_h5ad_dimensions(f_path_csc)

expect_equal(
  current = h5_meta_csr$dims,
  target = expected_dims,
  info = paste("h5ad ingestion - correct dimensions - CSR format")
)

expect_equal(
  current = h5_meta_csc$dims,
  target = expected_dims,
  info = paste("h5ad ingestion - correct dimensions - CSC format")
)

expect_equal(
  current = "CSR",
  target = h5_meta_csr$type,
  info = paste("h5ad ingestion - correct format - CSR format")
)

expect_equal(
  current = "CSC",
  target = h5_meta_csc$type,
  info = paste("h5ad ingestion - correct format - CSC format")
)

## loading ---------------------------------------------------------------------

### direct loading via helpers -------------------------------------------------

# test the loading helpers in rust

direct_load_csr <- rs_h5ad_data(
  f_path = f_path_csr,
  cs_type = h5_meta_csr$type,
  nrows = h5_meta_csr$dims[1],
  ncols = h5_meta_csr$dims[2],
  cell_quality = params_sc_min_quality(
    min_unique_genes = 0L,
    min_lib_size = 0L,
    min_cells = 0L,
    target_size = 1e5
  ),
  verbose = FALSE
)

direct_load_csc <- rs_h5ad_data(
  f_path = f_path_csc,
  cs_type = h5_meta_csc$type,
  nrows = h5_meta_csc$dims[1],
  ncols = h5_meta_csc$dims[2],
  cell_quality = params_sc_min_quality(
    min_unique_genes = 0L,
    min_lib_size = 0L,
    min_cells = 0L,
    target_size = 1e5
  ),
  verbose = FALSE
)

expect_true(
  current = all(direct_load_csr$indptr == direct_load_csc$indptr),
  info = paste("h5ad ingestion - flipping from CSC to CSR works - indptr")
)

expect_true(
  current = all(direct_load_csr$indices == direct_load_csc$indices),
  info = paste("h5ad ingestion - flipping from CSC to CSR works - indices")
)

expect_true(
  current = all(direct_load_csr$data == direct_load_csc$data),
  info = paste("h5ad ingestion - flipping from CSC to CSR works - data")
)

expect_true(
  current = all(direct_load_csr$indptr == single_cell_test_data$counts@p),
  info = paste("h5ad ingestion - indptr are preserved")
)

expect_true(
  current = all(direct_load_csr$indices == single_cell_test_data$counts@j),
  info = paste("h5ad ingestion - indices are preserved")
)

expect_true(
  current = all(direct_load_csr$data == single_cell_test_data$counts@x),
  info = paste("h5ad ingestion - data are preserved")
)

### qc params work -------------------------------------------------------------

#### expected data -------------------------------------------------------------

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
  info = "h5 - sensible amount of genes pass"
)

expect_true(
  current = length(cells_pass) > 800 & length(cells_pass) != 1000,
  info = "h5 - sensible amount of cells pass"
)

counts_filtered <- single_cell_test_data$counts[cells_pass, genes_pass]

counts_filtered_csc <- counts_csc[cells_pass, genes_pass]

obs_filtered <- single_cell_test_data$obs[cells_pass, ]

vars_filtered <- single_cell_test_data$var[genes_pass, ]

#### rust logic ----------------------------------------------------------------

sc_qc_param = params_sc_min_quality(
  min_unique_genes = min_genes_exp,
  min_lib_size = min_lib_size,
  min_cells = min_cells_exp
)

#### direct writing ------------------------------------------------------------

# test the underlying rust directly
sc_object <- SingleCells(dir_data = test_temp_dir)

rust_con <- get_sc_rust_ptr(sc_object)

file_res <- rust_con$h5_to_file(
  cs_type = h5_meta_csr$type,
  h5_path = path.expand(f_path_csr),
  no_cells = h5_meta_csr$dims["obs"],
  no_genes = h5_meta_csr$dims["var"],
  qc_params = sc_qc_param,
  verbose = FALSE
)

expect_equivalent(
  current = file_res$cell_indices + 1,
  target = cells_pass,
  info = paste("h5ad to binary - correct cells being kept")
)

expect_equivalent(
  current = file_res$gene_indices + 1,
  target = genes_pass,
  info = paste("h5ad to binary - correct genes being kept")
)

file_res <- rust_con$h5_to_file(
  cs_type = h5_meta_csc$type,
  h5_path = path.expand(f_path_csc),
  no_cells = h5_meta_csc$dims["obs"],
  no_genes = h5_meta_csc$dims["var"],
  qc_params = sc_qc_param,
  verbose = FALSE
)

expect_equivalent(
  current = file_res$cell_indices + 1,
  target = cells_pass,
  info = paste("h5ad to binary - correct cells being kept (with reading CSC)")
)

expect_equivalent(
  current = file_res$gene_indices + 1,
  target = genes_pass,
  info = paste("h5ad to binary - correct genes being kept (with reading CSC)")
)

#### streaming engine ----------------------------------------------------------

file_res <- rust_con$h5_to_file_streaming(
  cs_type = h5_meta_csr$type,
  h5_path = path.expand(f_path_csr),
  no_cells = h5_meta_csr$dims["obs"],
  no_genes = h5_meta_csr$dims["var"],
  qc_params = sc_qc_param,
  verbose = FALSE
)

expect_equivalent(
  current = file_res$cell_indices + 1,
  target = cells_pass,
  info = paste("h5ad to binary - correct cells being kept")
)

expect_equivalent(
  current = file_res$gene_indices + 1,
  target = genes_pass,
  info = paste("h5ad to binary - correct genes being kept")
)

file_res <- rust_con$h5_to_file_streaming(
  cs_type = h5_meta_csc$type,
  h5_path = path.expand(f_path_csc),
  no_cells = h5_meta_csc$dims["obs"],
  no_genes = h5_meta_csc$dims["var"],
  qc_params = sc_qc_param,
  verbose = FALSE
)

expect_equivalent(
  current = file_res$cell_indices + 1,
  target = cells_pass,
  info = paste("h5ad to binary - correct cells being kept (with reading CSC)")
)

expect_equivalent(
  current = file_res$gene_indices + 1,
  target = genes_pass,
  info = paste("h5ad to binary - correct genes being kept (with reading CSC)")
)

### duckdb part ----------------------------------------------------------------

duckdb_con <- get_sc_duckdb(sc_object)

# test csr reading in

duckdb_con$populate_obs_from_h5(
  h5_path = path.expand(f_path_csr),
  filter = as.integer(file_res$cell_indices + 1)
)

duckdb_con$populate_vars_from_h5(
  h5_path = path.expand(f_path_csr),
  filter = as.integer(file_res$gene_indices + 1)
)

obs_filtered_db <- duckdb_con$get_obs_table()

vars_filtered_db <- duckdb_con$get_vars_table()

expect_equivalent(
  current = obs_filtered_db$cell_id,
  target = obs_filtered$cell_id,
  info = "obs table from h5ad csr - correct cells kept"
)

expect_equivalent(
  current = vars_filtered_db$gene_id,
  target = vars_filtered$gene_id,
  info = "obs table from h5ad csr - correct genes kept"
)

# test csc

duckdb_con$populate_obs_from_h5(
  h5_path = path.expand(f_path_csc),
  filter = as.integer(file_res$cell_indices + 1)
)

duckdb_con$populate_vars_from_h5(
  h5_path = path.expand(f_path_csc),
  filter = as.integer(file_res$gene_indices + 1)
)

obs_filtered_db <- duckdb_con$get_obs_table()

vars_filtered_db <- duckdb_con$get_vars_table()

expect_equivalent(
  current = obs_filtered_db$cell_id,
  target = obs_filtered$cell_id,
  info = "obs table from h5ad csr - correct cells kept"
)

expect_equivalent(
  current = vars_filtered_db$gene_id,
  target = vars_filtered$gene_id,
  info = "obs table from h5ad csr - correct genes kept"
)

#### gene /cell maps -----------------------------------------------------------

cell_map <- duckdb_con$get_obs_index_map()

gene_map <- duckdb_con$get_var_index_map()


expect_true(
  current = checkmate::checkInteger(
    cell_map,
    names = "named",
    len = nrow(obs_filtered)
  ),
  info = "cell map as expected"
)

expect_true(
  current = checkmate::checkInteger(
    gene_map,
    names = "named",
    len = nrow(vars_filtered)
  ),
  info = "gene map as expected"
)

## direct object load ----------------------------------------------------------

sc_object <- SingleCells(dir_data = test_temp_dir)

sc_object <- load_h5ad(
  object = sc_object,
  h5_path = path.expand(f_path_csr),
  sc_qc_param = sc_qc_param,
  .verbose = FALSE
)

### test obs -------------------------------------------------------------------

obs_object <- sc_object[[]]

expect_equal(
  current = obs_object$cell_id,
  target = obs_filtered$cell_id,
  info = "obs table from h5ad csr (from object) - correct cells kept"
)

expect_equal(
  current = obs_object$cell_grp,
  target = obs_filtered$cell_grp,
  info = "obs table from h5ad csr (from object) - cell group correct"
)

expect_equivalent(
  current = Matrix::rowSums(counts_filtered),
  target = obs_object$lib_size,
  info = "obs table from h5ad csr (from object) - library size correct"
)

expect_equivalent(
  current = Matrix::rowSums(counts_filtered != 0),
  target = obs_object$nnz,
  info = "obs table from h5ad csr (from object) - nnz correct"
)

### test vars ------------------------------------------------------------------

vars_object <- get_sc_var(sc_object)

expect_equal(
  current = vars_object$gene_id,
  target = vars_filtered$gene_id,
  info = "var table from h5ad csr (from object) - correct genes kept"
)

expect_true(
  current = checkmate::qtest(vars_object$no_cells_exp, "I+"),
  info = "gene NNZ loaded via h5ad direct load"
)

### test counts ----------------------------------------------------------------

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

### test getters ---------------------------------------------------------------

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

## streaming h5ad --------------------------------------------------------------

sc_object <- SingleCells(dir_data = test_temp_dir)

sc_object <- stream_h5ad(
  object = sc_object,
  h5_path = path.expand(f_path_csr),
  sc_qc_param = sc_qc_param,
  .verbose = FALSE
)

### test obs -------------------------------------------------------------------

obs_object <- sc_object[[]]

expect_equal(
  current = obs_object$cell_id,
  target = obs_filtered$cell_id,
  info = "obs table from h5ad csr (from object) - correct cells kept"
)

expect_equal(
  current = obs_object$cell_grp,
  target = obs_filtered$cell_grp,
  info = "obs table from h5ad csr (from object) - cell group correct"
)

expect_equivalent(
  current = Matrix::rowSums(counts_filtered),
  target = obs_object$lib_size,
  info = "obs table from h5ad csr (from object) - library size correct"
)

expect_equivalent(
  current = Matrix::rowSums(counts_filtered != 0),
  target = obs_object$nnz,
  info = "obs table from h5ad csr (from object) - nnz correct"
)

### test vars ------------------------------------------------------------------

vars_object <- get_sc_var(sc_object)

expect_equal(
  current = vars_object$gene_id,
  target = vars_filtered$gene_id,
  info = "var table from h5ad csr (from object) - correct genes kept"
)

expect_true(
  current = checkmate::qtest(vars_object$no_cells_exp, "I+"),
  info = "gene NNZ loaded via h5ad streaming load"
)

### test counts ----------------------------------------------------------------

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

### test getters ---------------------------------------------------------------

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

## multi file read -------------------------------------------------------------

### write to files to disk -----------------------------------------------------

#### file 1 --------------------------------------------------------------------

single_cell_test_data_1 <- generate_single_cell_test_data(seed = 1L)

f_path_csr_1 = file.path(test_temp_dir, "csr_test_1.h5ad")

write_h5ad_sc(
  f_path = f_path_csr_1,
  counts = single_cell_test_data_1$counts,
  obs = single_cell_test_data_1$obs,
  var = single_cell_test_data_1$var,
  .verbose = FALSE
)

#### file 2 --------------------------------------------------------------------

single_cell_test_data_2 <- generate_single_cell_test_data(seed = 2L)

f_path_csr_2 = file.path(test_temp_dir, "csr_test_2.h5ad")

write_h5ad_sc(
  f_path = f_path_csr_2,
  counts = single_cell_test_data_2$counts,
  obs = single_cell_test_data_2$obs,
  var = single_cell_test_data_2$var,
  .verbose = FALSE
)

### test file scanning ---------------------------------------------------------

h5ad_files <- list.files(test_temp_dir)

h5ad_files <- h5ad_files[
  grepl("test_1", h5ad_files) | grepl("test_2", h5ad_files)
]

h5ad_files_final <- file.path(test_temp_dir, h5ad_files)
names(h5ad_files_final) <- c("exp1", "exp2")

h5_tasks <- prescan_h5ad_files(h5_paths = h5ad_files_final)

expect_true(
  current = checkmate::testList(h5_tasks),
  info = "h5_task_list correct type"
)

expect_true(
  current = all(
    c("universe", "universe_size", "file_tasks") %in%
      names(h5_tasks)
  ),
  info = "h5_task_list expected names"
)

### test rust part directly ----------------------------------------------------

# test the underlying rust directly
sc_object <- SingleCells(dir_data = test_temp_dir)

rust_con <- get_sc_rust_ptr(sc_object)

file_res <- rust_con$multi_h5_to_file(
  file_tasks = h5_tasks$file_tasks,
  universe_size = as.integer(h5_tasks$universe_size),
  qc_params = sc_qc_param,
  verbose = FALSE
)

expect_true(
  current = checkmate::testList(file_res),
  info = "h5 rust file rult is a list"
)

expect_true(
  current = checkmate::testNames(
    names(file_res),
    must.include = c(
      "global_gene_indices",
      "total_cells",
      "total_genes",
      "per_file"
    )
  ),
  info = "h5 rust file result has expected structure"
)

counts <- rust_con$return_full_mat(
  assay = "raw",
  cell_based = FALSE,
  verbose = FALSE
)

expect_true(
  current = counts$no_cells > 1000,
  info = "more than 1000 cells were written to file"
)

### test duckdb part -----------------------------------------------------------

#### obs -----------------------------------------------------------------------

duckdb_con <- get_sc_duckdb(sc_object)

per_file_info <- lapply(file_res$per_file, function(f) {
  list(
    h5_path = h5_tasks$file_tasks[[f$exp_id]]$h5_path,
    exp_id = f$exp_id,
    cell_filter = as.integer(f$cell_indices + 1L)
  )
})

duckdb_con$populate_obs_from_multi_h5(
  per_file_info = per_file_info,
  cell_id_col = NULL
)

obs_direct <- duckdb_con$get_obs_table()

expect_true(
  current = checkmate::testDataTable(obs_direct),
  info = "obs table correctly being written"
)

expect_true(
  current = nrow(obs_direct) == counts$no_cells,
  info = "counts and obs match in multi-file import"
)

#### vars ----------------------------------------------------------------------

final_gene_names <- h5_tasks$universe[file_res$global_gene_indices + 1L]

duckdb_con$populate_vars_from_h5_reordered(
  h5_path = h5_tasks$file_tasks[[1L]]$h5_path,
  final_gene_names = final_gene_names
)

var_direct <- duckdb_con$get_vars_table()

expect_true(
  current = checkmate::testDataTable(var_direct),
  info = "var table correctly being written"
)

expect_true(
  current = nrow(var_direct) == counts$no_genes,
  info = "counts and vars match in multi-file import"
)

### main method ----------------------------------------------------------------

sc_object <- load_multi_h5ad(
  object = sc_object,
  prescan_result = h5_tasks,
  sc_qc_param = sc_qc_param,
  .verbose = FALSE
)

expect_equal(
  current = sc_object[[colnames(obs_direct)]],
  target = obs_direct,
  info = "main multi-read method returns an obs table"
)

expect_equal(
  current = dim(sc_object[]),
  target = c(counts$no_cells, counts$no_genes),
  info = "main multi-read method returns counts"
)

expect_equal(
  current = get_sc_var(sc_object)[, colnames(var_direct), with = FALSE],
  target = var_direct,
  info = "main multi-read method returns a var table"
)

# clean up ---------------------------------------------------------------------

on.exit(unlink(test_temp_dir, recursive = TRUE, force = TRUE), add = TRUE)
