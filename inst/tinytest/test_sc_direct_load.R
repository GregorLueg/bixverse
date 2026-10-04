# generate synthetic data ------------------------------------------------------

source("helper_sc.R", local = TRUE)

library(magrittr)

test_temp_dir <- sc_test_dir("direct_load")

## params ----------------------------------------------------------------------

# thresholds
min_lib_size <- 300L
min_genes_exp <- 45L
min_cells_exp <- 500L

## synthetic data --------------------------------------------------------------

# test the load from disk

single_cell_test_data <- sc_test_fixture(
  min_lib_size = min_lib_size,
  min_genes_exp = min_genes_exp,
  min_cells_exp = min_cells_exp
)

## generate the object ---------------------------------------------------------

sc_qc_param <- sc_test_qc_params(single_cell_test_data, target_size = 1000)

sc_object <- sc_test_object(
  test_temp_dir,
  single_cell_test_data,
  sc_qc_param = sc_qc_param
)

# do a filtering on the obs column
sc_object <- set_cells_to_keep(sc_object, unlist(sc_object[["cell_id"]][1:500]))

# check if the NNZ per gene was added

expect_true(
  current = checkmate::qtest(get_sc_var(sc_object)[["no_cells_exp"]], "I+"),
  info = "gene NNZ added by R direct load"
)

# remove it...
rm(sc_object)

# tests ------------------------------------------------------------------------

## load from disk --------------------------------------------------------------

sc_object <- SingleCells(dir_data = test_temp_dir)

sc_object <- load_existing(sc_object, .verbose = FALSE)

## getter checks ---------------------------------------------------------------

expect_true(
  current = checkmate::qtest(get_cell_names(sc_object), "S+"),
  info = "loading from disk directly - cell names"
)

expect_true(
  current = checkmate::qtest(get_gene_names(sc_object), "S+"),
  info = "loading from disk directly - gene names"
)

obs_dt <- sc_object[[]]
var_dt <- get_sc_var(sc_object)

expect_true(
  current = checkmate::testDataTable(obs_dt),
  info = "loading from disk directly - obs table"
)

expect_true(
  current = checkmate::testDataTable(var_dt),
  info = "loading from disk directly - var table"
)

expect_true(
  current = nrow(obs_dt) == 500L,
  info = "obs_table filtering works"
)

expect_true(
  current = nrow(get_sc_obs(sc_object)) == sc_object@dims[1],
  info = "obs_table full table is still accessible"
)

expect_true(
  current = length(get_cells_to_keep(sc_object)) == 500L,
  info = "cell_to_keep filtering also worked"
)

## count getters ---------------------------------------------------------------

cell_counts <- sc_object[]

gene_counts <- sc_object[,, return_format = "gene"]

expect_true(
  current = checkmate::testClass(cell_counts, "dgRMatrix"),
  info = "loading from disk directly - counts in cell friendly format"
)

expect_true(
  current = checkmate::testClass(gene_counts, "dgCMatrix"),
  info = "loading from disk directly - counts in gene friendly format"
)

expect_true(
  current = dim(cell_counts)[2] == nrow(var_dt),
  info = "loading from disk directly - expected dimensions of counts and var"
)

expect_equivalent(
  current = dim(gene_counts),
  target = dim(cell_counts),
  info = "same dimensions for the counts"
)

## save file with memory cached data -------------------------------------------

# add data to the object
sc_object <- find_hvg_sc(sc_object, hvg_no = 30L, .verbose = FALSE)

sc_object <- calculate_pca_sc(
  sc_object,
  no_pcs = 5L,
  .verbose = FALSE
)

hvg_genes_initial <- get_hvg(sc_object)
pca_factors_initial <- get_pca_factors(sc_object)

### saving to disk -------------------------------------------------------------

save_sc_exp_to_disk(sc_object)

save_sc_exp_to_disk(sc_object, type = "rds")

expect_true(
  current = "memory.rds" %in% list.files(path = test_temp_dir),
  info = "RDS saving works"
)

expect_true(
  current = "memory.qs2" %in% list.files(path = test_temp_dir),
  info = "qs2 saving works"
)

### qs2 ------------------------------------------------------------------------

rm(sc_object)

sc_object <- SingleCells(dir_data = test_temp_dir)

expect_message(current = load_existing(sc_object), info = "message working")

sc_object <- load_existing(sc_object, .verbose = FALSE)

expect_equal(
  current = get_pca_factors(sc_object),
  target = pca_factors_initial,
  info = "PCA loaded in correctly - qs2"
)

expect_equal(
  current = get_hvg(sc_object),
  target = hvg_genes_initial,
  info = "HVGs loaded in correctly - qs2"
)

### rds ------------------------------------------------------------------------

rm(sc_object)

sc_object <- SingleCells(dir_data = test_temp_dir)

# will force the function to load from rds
removed <- file.remove(file.path(test_temp_dir, "memory.qs2"))

sc_object <- load_existing(sc_object, .verbose = FALSE)

expect_equal(
  current = get_pca_factors(sc_object),
  target = pca_factors_initial,
  info = "PCA loaded in correctly - RDS"
)

expect_equal(
  current = get_hvg(sc_object),
  target = hvg_genes_initial,
  info = "HVGs loaded in correctly - RDS"
)

### archive --------------------------------------------------------------------

raw_before <- sc_object[]
norm_before <- sc_object[,, assay = "norm"]
genes_before <- sc_object[,, return_format = "gene"]

archive_stats <- archive_sc_exp(sc_object, level = 3L, .verbose = FALSE)

expect_true(
  current = file.exists(file.path(test_temp_dir, "counts.bxa")) &&
    !any(file.exists(file.path(
      test_temp_dir,
      c("counts_cells.bin", "counts_genes.bin")
    ))),
  info = "archive written and binaries removed"
)

expect_equal(
  current = archive_stats$n_norm_stored,
  target = 0L,
  info = "all norms recompute from the raw counts"
)

rm(sc_object)

sc_object <- load_existing(
  SingleCells(dir_data = test_temp_dir),
  .verbose = FALSE
)

expect_equal(
  current = sc_object[],
  target = raw_before,
  info = "raw counts survive the archive round trip"
)

expect_equal(
  current = sc_object[,, assay = "norm"],
  target = norm_before,
  info = "norm counts survive the archive round trip"
)

expect_equal(
  current = sc_object[,, return_format = "gene"],
  target = genes_before,
  info = "gene file is rebuilt on restore"
)

# clean up ---------------------------------------------------------------------

sc_test_cleanup(test_temp_dir)
