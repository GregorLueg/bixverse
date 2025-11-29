# differential abundance tests -------------------------------------------------

## test parameters -------------------------------------------------------------

# thresholds
min_lib_size <- 300L
min_genes_exp <- 45L
min_cells_exp <- 500L
hvg_to_keep <- 50L
no_pcs <- 20L

test_temp_dir <- file.path(
  tempdir(),
  paste0("test_", format(Sys.time(), "%Y%m%d_%H%M%S_"), sample(1000:9999, 1))
)
dir.create(test_temp_dir, recursive = TRUE)

## synthetic test data ---------------------------------------------------------

single_cell_test_data <- generate_single_cell_test_data(
  syn_data_params = params_sc_synthetic_data(
    n_samples = 6L,
    sample_bias = "even"
  )
)

# helper function to not have to regenerate the synthetic data...
sample_ids_uneven <- bixverse:::rs_sample_ids_for_cell_types(
  cell_type_indices = as.integer(factor(single_cell_test_data$obs$cell_grp)) -
    1L,
  n_samples = 6,
  sample_bias = "very_uneven",
  seed = 123L
)

sample_ids_uneven <- sprintf("sample_%i", sample_ids_uneven + 1L)

single_cell_test_data$obs[, sample_id_v2 := sample_ids_uneven]

design_df <- data.frame(
  grps = c(rep("ctr", 3), rep("trt", 3))
)
rownames(design_df) <- sprintf("sample_%i", 1:6)

# sensible amount of genes and cells passing the thresholds
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

## object gen ------------------------------------------------------------------

sc_object <- single_cell_exp(dir_data = test_temp_dir)

sc_object <- load_r_data(
  object = sc_object,
  counts = single_cell_test_data$counts,
  obs = single_cell_test_data$obs,
  var = single_cell_test_data$var,
  sc_qc_param = params_sc_min_quality(
    min_unique_genes = min_genes_exp,
    min_lib_size = min_lib_size,
    min_cells = min_cells_exp,
    target_size = 1000
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
  sc_object,
  no_pcs = no_pcs,
  .verbose = FALSE,
  randomised_svd = TRUE
)

# tests ------------------------------------------------------------------------

## warnings --------------------------------------------------------------------

expect_warning(
  current = get_miloR_abundances_sc(
    object = sc_object,
    sample_id_col = "sample_id",
    .verbose = FALSE
  ),
  info = "miloR warnings if there is no kNN generated"
)

sc_object <- find_neighbours_sc(sc_object, .verbose = FALSE)

expect_warning(
  current = get_miloR_abundances_sc(
    object = sc_object,
    sample_id_col = "sample_id",
    miloR_params = params_sc_miloR(k_refine = 100L),
    .verbose = FALSE
  ),
  info = "miloR complains if k_refine >= no embedding dimensions"
)

## neighbourhood generation ----------------------------------------------------

miloR_obj_index <- get_miloR_abundances_sc(
  object = sc_object,
  sample_id_col = "sample_id",
  miloR_params = params_sc_miloR(refinement_strategy = "index"),
  .verbose = FALSE
)

miloR_obj_approx <- get_miloR_abundances_sc(
  object = sc_object,
  sample_id_col = "sample_id",
  miloR_params = params_sc_miloR(refinement_strategy = "approximate"),
  .verbose = FALSE
)

miloR_obj_bruteforce <- get_miloR_abundances_sc(
  object = sc_object,
  sample_id_col = "sample_id",
  miloR_params = params_sc_miloR(refinement_strategy = "bruteforce"),
  .verbose = FALSE
)

expect_inherits(
  current = miloR_obj_index,
  class = "sc_miloR",
  info = "correct class returned index strategy"
)

expect_inherits(
  current = miloR_obj_approx,
  class = "sc_miloR",
  info = "correct class returned approx strategy"
)

expect_inherits(
  current = miloR_obj_bruteforce,
  class = "sc_miloR",
  info = "correct class returned brute force strategy"
)

index_cells_index <- as.character(get_index_cells(miloR_obj_index))
index_cells_approx <- as.character(get_index_cells(miloR_obj_approx))
index_cells_bruteforce <- as.character(get_index_cells(miloR_obj_bruteforce))

expect_true(
  current = rs_set_similarity(index_cells_index, index_cells_bruteforce, TRUE) >
    0.8,
  info = "high overlap in index cells for index method and brute force method"
)

# this will yield lower overlaps due to the constraints in the approx method
expect_true(
  current = rs_set_similarity(index_cells_index, index_cells_approx, TRUE) >
    0.6,
  info = "lower overlap in index cells for approx method and approx method"
)

expect_true(
  current = rs_set_similarity(
    index_cells_approx,
    index_cells_bruteforce,
    TRUE
  ) >
    0.6,
  info = "lower overlap in index cells for approx method and brute force method"
)

# spatial distances are sensible
expect_true(
  all(miloR_obj_index$spatial_dist > 0),
  info = "all spatial distances should be positive"
)
expect_true(
  all(is.finite(miloR_obj_index$spatial_dist)),
  info = "all spatial distances should be finite"
)
expect_true(
  median(miloR_obj_index$spatial_dist) > 0,
  info = "median spatial distance should be positive"
)

# expected dimensions in sample counts
expect_equal(
  current = ncol(miloR_obj_index$sample_counts),
  target = nrow(unique(sc_object[["sample_id"]])),
  info = "sample counts should have one column per sample"
)

# correct params stored
expect_equal(
  miloR_obj_index$params$refinement_strategy,
  target = "index",
  info = "correct parameters stored - index version"
)
expect_equal(
  miloR_obj_approx$params$refinement_strategy,
  target = "approximate",
  info = "correct parameters stored - approximate version"
)
expect_equal(
  miloR_obj_bruteforce$params$refinement_strategy,
  target = "bruteforce",
  info = "correct parameters stored - bruteforce version"
)

# check that the proportions are having
miloR_obj_index_small <- get_miloR_abundances_sc(
  object = sc_object,
  sample_id_col = "sample_id",
  miloR_params = params_sc_miloR(refinement_strategy = "index", prop = 0.1),
  .verbose = FALSE
)

expect_true(
  current = length(get_index_cells(miloR_obj_index_small)) <
    length(get_index_cells(miloR_obj_index)),
  info = "less index cells found with smaller proportion"
)

## neighbourhood testing -------------------------------------------------------

# this one SHOULD have significant ones
miloR_obj_index_v2 <- get_miloR_abundances_sc(
  object = sc_object,
  sample_id_col = "sample_id_v2",
  .verbose = FALSE
)

expect_equal(
  current = get_index_cells(miloR_obj_index_v2),
  target = get_index_cells(miloR_obj_index),
  info = "given a different sample id column there is no change"
)

miloR_obj_index <- test_nhoods(
  x = miloR_obj_index,
  design = ~grps,
  design_df = design_df
)

miloR_obj_index_v2 <- test_nhoods(
  x = miloR_obj_index_v2,
  design = ~grps,
  design_df = design_df
)

expect_true(
  current = checkmate::testDataTable(get_differential_abundance_res(
    miloR_obj_index
  )),
  info = "getter behaving and returning the right format"
)

expect_true(
  current = checkmate::testNames(
    names(get_differential_abundance_res(miloR_obj_index)),
    must.include = c(
      "Nhood",
      "logFC",
      "logCPM",
      "F",
      "PValue",
      "FDR",
      "SpatialFDR"
    )
  ),
  info = "expected columns on the results"
)

expect_true(
  current = !any(get_differential_abundance_res(miloR_obj_index)$FDR <= 0.05),
  info = "no significant differences found in the correct case"
)

expect_true(
  current = any(get_differential_abundance_res(miloR_obj_index_v2)$FDR <= 0.05),
  info = "significant differences found in the correct case"
)

miloR_obj_index <- test_nhoods(
  x = miloR_obj_index,
  design = ~grps,
  design_df = design_df,
  fdr_weighting = "none"
)

expect_true(
  current = all(is.na(
    get_differential_abundance_res(miloR_obj_index)$SpatialFDR
  )),
  info = "spatial FDR NOT calculated"
)

## neighbour information -------------------------------------------------------

# this one does not have the neighbourhood data, so should warn

expect_warning(
  current = add_nhoods_info(
    x = miloR_obj_index_small,
    cell_info = sc_object[[]]$cell_grp
  ),
  info = paste(
    "warning that addition of cell_info does not work if there is",
    "no neighbourhood information"
  )
)

miloR_obj_index <- add_nhoods_info(
  x = miloR_obj_index,
  cell_info = sc_object[[]]$cell_grp
)

expect_true(
  current = all(miloR_obj_index$nhoods_info$majority_prop > 0.9),
  info = "majority of the neighbourhoods are the same cell type"
)

on.exit(unlink(test_temp_dir, recursive = TRUE, force = TRUE), add = TRUE)
