# generate synthetic data ------------------------------------------------------

## params ----------------------------------------------------------------------

# thresholds
min_lib_size <- 300L
min_genes_exp <- 45L
min_cells_exp <- 500L

## synthetic data --------------------------------------------------------------

# test the load from disk

single_cell_test_data <- generate_single_cell_test_data()

f_path_csr = file.path(tempdir(), "csr_test.h5ad")

write_h5ad_sc(
  f_path = f_path_csr,
  counts = single_cell_test_data$counts,
  obs = single_cell_test_data$obs,
  single_cell_test_data$var,
  .verbose = FALSE
)

## generate the object ---------------------------------------------------------

sc_qc_param = params_sc_min_quality(
  min_unique_genes = min_genes_exp,
  min_lib_size = min_lib_size,
  min_cells = min_cells_exp,
  target_size = 1000
)

sc_object <- suppressWarnings(single_cell_exp(dir_data = tempdir()))

sc_object <- load_h5ad(
  object = sc_object,
  h5_path = path.expand(f_path_csr),
  sc_qc_param = sc_qc_param,
  .verbose = FALSE
)

# do a filtering on the obs column
sc_object <- set_cells_to_keep(sc_object, unlist(sc_object[["cell_id"]][1:500]))

# remove it...
rm(sc_object)

# tests ------------------------------------------------------------------------

## load from disk --------------------------------------------------------------

sc_object <- single_cell_exp(dir_data = tempdir())

sc_object <- load_existing(sc_object)

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
