# generate synthetic data ------------------------------------------------------

test_temp_dir <- file.path(
  tempdir(),
  paste0("test_", format(Sys.time(), "%Y%m%d_%H%M%S_"), sample(1000:9999, 1))
)
dir.create(test_temp_dir, recursive = TRUE)

## params ----------------------------------------------------------------------

# thresholds
min_lib_size <- 300L
min_genes_exp <- 45L
min_cells_exp <- 500L

## synthetic data --------------------------------------------------------------

# test the load from disk

single_cell_test_data <- generate_single_cell_test_data()

## generate the object ---------------------------------------------------------

sc_qc_param = params_sc_min_quality(
  min_unique_genes = min_genes_exp,
  min_lib_size = min_lib_size,
  min_cells = min_cells_exp,
  target_size = 1000
)

sc_object <- single_cell_exp(dir_data = test_temp_dir)

sc_object <- load_r_data(
  object = sc_object,
  counts = single_cell_test_data$counts,
  obs = single_cell_test_data$obs,
  var = single_cell_test_data$var,
  sc_qc_param = sc_qc_param,
  .verbose = FALSE
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

sc_object <- single_cell_exp(dir_data = test_temp_dir)

sc_object <- suppressMessages(load_existing(sc_object))

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

sc_object <- calculate_pca_sc(sc_object, no_pcs = 5L, .verbose = FALSE)

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

sc_object <- single_cell_exp(dir_data = test_temp_dir)

expect_message(current = load_existing(sc_object), info = "message working")

sc_object <- suppressMessages(load_existing(sc_object))

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

sc_object <- single_cell_exp(dir_data = test_temp_dir)

# will force the function to load from rds
removed <- file.remove(file.path(test_temp_dir, "memory.qs2"))

sc_object <- suppressMessages(load_existing(sc_object))

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

on.exit(unlink(test_temp_dir, recursive = TRUE, force = TRUE), add = TRUE)
