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
  neighbours_params = params_sc_neighbours(knn = list(k = 15L)),
  .verbose = FALSE
)

## meta cells ------------------------------------------------------------------

### hdwgcna --------------------------------------------------------------------

hdwgcna <- get_meta_cells_sc(
  sc_object,
  sc_meta_cell_params = params_sc_metacells(
    target_no_metacells = 50L
  ),
  .verbose = FALSE
)

hdwgcna <- calc_meta_cell_purity(
  hdwgcna,
  original_cell_type = unlist(sc_object[["cell_grp"]])
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


expect_true(
  current = mean(hdwgcna[[]]$mc_purity) > 0.75,
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

#### subsetted version ---------------------------------------------------------

# to imitate a scenario in which someone would like to generate meta cells
# in specific cell types; will remove cell_type_3
cells_to_use <- sc_object[[]][cell_grp != "cell_type_3", cell_id]

hdwgcna_small <- get_meta_cells_sc(
  sc_object,
  sc_meta_cell_params = params_sc_metacells(
    target_no_metacells = 50L,
    max_shared = 10L,
    knn = list(k = 10L)
  ),
  cells_to_use = cells_to_use,
  .verbose = FALSE
)

hdwgcna_small <- calc_meta_cell_purity(
  hdwgcna_small,
  original_cell_type = unlist(sc_object[["cell_grp"]])
)

original_cell_type = unlist(sc_object[["cell_grp"]])

right_cell_types <- purrr::map_lgl(
  hdwgcna_small[[]]$original_cell_idx,
  function(idx) {
    types <- original_cell_type[idx]
    any(types != "cell_type_3")
  }
)

expect_true(
  current = mean(hdwgcna_small[[]]$mc_purity) > 0.75,
  info = paste(
    "hgwgnca - similar cell types are being pulled together;",
    "subsetted version"
  )
)

expect_true(
  current = all(right_cell_types),
  info = "no unexpected cell types in subsetted version - hdwgcna"
)


### seacells -------------------------------------------------------------------

seacells <- get_seacells_sc(
  sc_object,
  seacell_params = params_sc_seacells(
    n_sea_cells = 50L,
    min_iter = 5L,
    knn = list(k = 10L)
  ),
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

seacells <- calc_meta_cell_purity(
  seacells,
  original_cell_type = unlist(sc_object[["cell_grp"]])
)

# much better than hdwgcna!
expect_true(
  current = mean(seacells[[]]$mc_purity) > 0.85,
  info = "seacell - similar cell types are being pulled together"
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

#### smaller subset ------------------------------------------------------------

seacells_small <- get_seacells_sc(
  sc_object,
  seacell_params = params_sc_seacells(
    n_sea_cells = 50L,
    min_iter = 5L,
    convergence_epsilon = 0.001,
    knn = list(k = 10L)
  ),
  cells_to_use = cells_to_use,
  .verbose = FALSE
)

seacells_small <- calc_meta_cell_purity(
  seacells_small,
  original_cell_type = unlist(sc_object[["cell_grp"]])
)

right_cell_types <- purrr::map_lgl(
  seacells_small[[]]$original_cell_idx,
  function(idx) {
    types <- original_cell_type[idx]
    any(types != "cell_type_3")
  }
)

expect_true(
  current = mean(seacells_small[[]]$mc_purity) > 0.75,
  info = paste(
    "hgwgnca - similar cell types are being pulled together;",
    "subsetted version"
  )
)

expect_true(
  current = all(right_cell_types),
  info = "no unexpected cell types in subsetted version - seacell"
)

### supercells -----------------------------------------------------------------

graining_factor <- 20

supercells <- get_supercells_sc(
  sc_object,
  sc_supercell_params = params_sc_supercell(
    graining_factor = graining_factor,
    knn = list(k = 10L)
  ),
  regenerate_knn = TRUE,
  .verbose = FALSE
)

expected_meta_cells <- ceiling(sc_object@dims[1] / graining_factor)

expect_true(
  current = checkmate::testDataTable(
    supercells[[]],
    nrows = expected_meta_cells
  ),
  info = "supercells obs correct - correct dimensions and type"
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
  info = "supercells obs correct - correct columns"
)

supercells <- calc_meta_cell_purity(
  supercells,
  original_cell_type = unlist(sc_object[["cell_grp"]])
)

expect_true(
  current = mean(supercells[[]]$mc_purity) > 0.85,
  info = "supercell - similar cell types are being pulled together"
)

expect_true(
  current = checkmate::testDataTable(
    get_sc_var(supercells)
  ),
  info = "supercell meta cells var correct"
)

expect_equal(
  current = dim(supercells[]),
  target = c(expected_meta_cells, 81),
  info = "supercell - correct return dimensions for the raw counts"
)

expect_true(
  current = checkmate::testClass(supercells[], "dgRMatrix"),
  info = "supercell - correct compressed sparse matrix type"
)

#### smaller subset ------------------------------------------------------------

supercell_small <- get_supercells_sc(
  sc_object,
  sc_supercell_params = params_sc_supercell(
    graining_factor = graining_factor,
    knn = list(k = 10L)
  ),
  cells_to_use = cells_to_use,
  .verbose = FALSE
)

supercell_small <- calc_meta_cell_purity(
  supercell_small,
  original_cell_type = unlist(sc_object[["cell_grp"]])
)

right_cell_types <- purrr::map_lgl(
  supercell_small[[]]$original_cell_idx,
  function(idx) {
    types <- original_cell_type[idx]
    any(types != "cell_type_3")
  }
)

expect_true(
  current = mean(supercell_small[[]]$mc_purity) > 0.75,
  info = paste(
    "supercell - similar cell types are being pulled together;",
    "subsetted version"
  )
)

expect_true(
  current = all(right_cell_types),
  info = "no unexpected cell types in subsetted version - supercell"
)

# clean up ---------------------------------------------------------------------

on.exit(unlink(test_temp_dir, recursive = TRUE, force = TRUE), add = TRUE)
