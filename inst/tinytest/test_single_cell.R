# h5 io ------------------------------------------------------------------------

## synthetic data --------------------------------------------------------------

single_cell_test_data <- generate_single_cell_test_data()

# CSR version
f_path_csr = file.path(tempdir(), "csr_test.h5ad")

write_h5ad_sc(
  f_path = f_path_csr,
  counts = single_cell_test_data$counts,
  obs = single_cell_test_data$obs,
  single_cell_test_data$var,
  .verbose = FALSE
)

# CSC version
counts_csc <- as(single_cell_test_data$counts, "CsparseMatrix")

f_path_csc = file.path(tempdir(), "csc_test.h5ad")

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

h5_meta_csr <- get_h5ad_dimensions(f_path_csr)

h5_meta_csc <- get_h5ad_dimensions(f_path_csc)

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
    target_size = 1000
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
    target_size = 1000
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

# these are a bit random, but at least I can check indices with this
min_lib_size <- 150L
min_genes_exp <- 10L
min_cells_exp <- 15L

genes_pass <- which(
  Matrix::colSums(single_cell_test_data$counts != 0) >= min_cells_exp
)

cells_pass <- which(
  (Matrix::rowSums(single_cell_test_data$counts[, genes_pass]) >=
    min_lib_size) &
    (Matrix::rowSums(single_cell_test_data$counts[, genes_pass] != 0) >=
      min_genes_exp)
)

sc_qc_param = params_sc_min_quality(
  min_unique_genes = min_genes_exp,
  min_lib_size = min_lib_size,
  min_cells = min_cells_exp
)

# test the underlying rust directly
sc_object <- suppressWarnings(single_cell_exp(dir_data = tempdir()))

rust_con <- get_sc_rust_ptr(sc_object)

file_res <- rust_con$h5_to_file(
  cs_type = h5_meta_csr$type,
  h5_path = path.expand(f_path_csr),
  no_cells = h5_meta_csr$dims["obs"],
  no_genes = h5_meta_csr$dims["var"],
  qc_params = sc_qc_param,
  verbose = TRUE
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
  verbose = TRUE
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
