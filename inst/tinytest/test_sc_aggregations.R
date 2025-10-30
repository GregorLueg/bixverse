# sc aggregations --------------------------------------------------------------

library(magrittr)

test_temp_dir <- file.path(
  tempdir(),
  paste0("test_", format(Sys.time(), "%Y%m%d_%H%M%S_"), sample(1000:9999, 1))
)

dir.create(test_temp_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot("Test directory does not exist" = dir.exists(test_temp_dir))

## testing parameters ----------------------------------------------------------

# thresholds
min_lib_size <- 300L
min_genes_exp <- 45L
min_cells_exp <- 500L
# hvg
hvg_to_keep <- 30L
# pca
no_pcs <- 10L

## synthetic test data ---------------------------------------------------------

single_cell_test_data <- generate_single_cell_test_data()

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
  info = "sc processing - sensible amount of genes pass"
)

expect_true(
  current = length(cells_pass) > 800 & length(cells_pass) != 1000,
  info = "sc processing - sensible amount of cells pass"
)

sc_qc_param <- params_sc_min_quality(
  min_unique_genes = min_genes_exp,
  min_lib_size = min_lib_size,
  min_cells = min_cells_exp,
  target_size = 1000
)

## underlying class ------------------------------------------------------------

sc_object <- single_cell_exp(dir_data = test_temp_dir)

sc_object <- # keep all cells for the sake of this
  sc_object <- load_r_data(
    object = sc_object,
    counts = single_cell_test_data$counts,
    obs = single_cell_test_data$obs,
    var = single_cell_test_data$var,
    sc_qc_param = params_sc_min_quality(
      min_unique_genes = min_genes_exp,
      min_lib_size = min_lib_size,
      min_cells = min_cells_exp
    ),
    streaming = FALSE,
    .verbose = FALSE
  )

sc_object <- find_hvg_sc(
  object = sc_object,
  hvg_no = hvg_to_keep,
  .verbose = FALSE
)

sc_object <- calculate_pca_sc(
  object = sc_object,
  no_pcs = no_pcs,
  randomised_svd = TRUE,
  .verbose = FALSE
)

sc_object <- find_neighbours_sc(
  sc_object,
  neighbours_params = params_sc_neighbours(knn_algorithm = "annoy"),
  .verbose = FALSE
)

## meta cells ------------------------------------------------------------------

cell_idx_to_type <- sc_object[[c("cell_idx", "cell_grp")]] %$%
  setNames(cell_grp, cell_idx)

### hdwgcna --------------------------------------------------------------------

hdwgcna <- get_meta_cells_sc(
  sc_object,
  sc_meta_cell_params = params_sc_metacells(target_no_metacells = 50L, k = 5L),
  .verbose = FALSE
)

expect_true(
  current = checkmate::testDataTable(
    hdwgcna[[]],
    nrows = 50
  ),
  info = "hdwgcna meta cells obs correct - correct dimensions and type"
)

expect_true(
  current = checkmate::testNames(
    names(hdwgcna[[]]),
    must.include = c(
      "meta_cell_idx",
      "meta_cell_id",
      "no_originating_cells",
      "original_cell_idx"
    )
  ),
  info = "hdwgcna meta cells obs correct - correct columns"
)

grouped_cell_types_max <- sapply(hdwgcna[[]]$original_cell_idx, function(idx) {
  types <- cell_idx_to_type[as.character(idx)]
  max(table(types)) / length(types)
})

expect_true(
  current = mean(grouped_cell_types_max) > 0.75,
  info = "similar cell types are being pulled together"
)

expect_true(
  current = checkmate::testDataTable(
    get_sc_var(hdwgcna)
  ),
  info = "hdwgcna meta cells var correct"
)

expect_equal(
  current = dim(hdwgcna[]),
  target = c(50, 81),
  info = "hdwgcna meta cells - correct return dimensions for the raw counts"
)

expect_true(
  current = checkmate::testClass(hdwgcna[], "dgRMatrix"),
  info = "hdwgcna meta cells - correct compressed sparse matrix type"
)

### seacells -------------------------------------------------------------------

seacells <- get_seacells_sc(
  sc_object,
  seacell_params = params_sc_seacells(n_sea_cells = 50L, k = 5L, min_iter = 5L),
  .verbose = FALSE
)

expect_true(
  current = checkmate::testDataTable(
    seacells[[]],
    nrows = 50
  ),
  info = "seacells obs correct - correct dimensions and type"
)

expect_true(
  current = checkmate::testNames(
    names(seacells[[]]),
    must.include = c(
      "meta_cell_idx",
      "meta_cell_id",
      "no_originating_cells",
      "original_cell_idx"
    )
  ),
  info = "seacells obs correct - correct columns"
)

grouped_cell_types_max <- sapply(seacells[[]]$original_cell_idx, function(idx) {
  types <- cell_idx_to_type[as.character(idx)]
  max(table(types)) / length(types)
})

# much better than hdwgcna!
expect_true(
  current = mean(grouped_cell_types_max) > 0.85,
  info = "similar cell types are being pulled together"
)

expect_true(
  current = checkmate::testDataTable(
    get_sc_var(seacells)
  ),
  info = "seacells meta cells var correct"
)

expect_equal(
  current = dim(seacells[]),
  target = c(50, 81),
  info = "seacells - correct return dimensions for the raw counts"
)

expect_true(
  current = checkmate::testClass(seacells[], "dgRMatrix"),
  info = "seacells - correct compressed sparse matrix type"
)

on.exit(unlink(test_temp_dir, recursive = TRUE, force = TRUE), add = TRUE)
