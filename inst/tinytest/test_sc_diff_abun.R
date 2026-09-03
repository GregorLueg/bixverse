# differential abundance tests -------------------------------------------------

source("helper_sc.R", local = TRUE)

library(magrittr)

test_temp_dir <- sc_test_dir("sc_diff_abund")

## test parameters -------------------------------------------------------------

# thresholds
min_lib_size <- 300L
min_genes_exp <- 45L
min_cells_exp <- 500L
hvg_to_keep <- 50L
no_pcs <- 20L
n_samples <- 6L

## synthetic test data ---------------------------------------------------------

single_cell_test_data <- sc_test_fixture(
  min_lib_size = min_lib_size,
  min_genes_exp = min_genes_exp,
  min_cells_exp = min_cells_exp,
  hvg_to_keep = hvg_to_keep,
  no_pcs = no_pcs,
  syn_data_params = params_sc_synthetic_data(
    n_samples = n_samples,
    sample_bias = "even"
  )
)

# helper function to not have to regenerate the synthetic data...
sample_ids_uneven <- bixverse:::rs_sample_ids_for_cell_types(
  cell_type_indices = as.integer(factor(single_cell_test_data$obs$cell_grp)) -
    1L,
  n_samples = n_samples,
  sample_bias = "very_uneven",
  seed = 123L
)

sample_ids_uneven <- sprintf("sample_%i", sample_ids_uneven + 1L)

single_cell_test_data$obs[, sample_id_v2 := sample_ids_uneven]

design_df <- data.frame(
  grps = c(rep("ctr", 3), rep("trt", 3))
)
rownames(design_df) <- sprintf("sample_%i", 1:6)

genes_pass <- single_cell_test_data$genes_pass
cells_pass <- single_cell_test_data$cells_pass

## object gen ------------------------------------------------------------------

sc_object <- sc_test_object(
  test_temp_dir,
  single_cell_test_data,
  sc_qc_param = sc_test_qc_params(single_cell_test_data, target_size = 1000)
)

sc_object <- find_hvg_sc(
  object = sc_object,
  hvg_no = hvg_to_keep,
  .verbose = FALSE
)

sc_object <- calculate_pca_sc(
  sc_object,
  no_pcs = no_pcs,
  .verbose = FALSE
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
  class = "miloR",
  info = "correct class returned index strategy"
)

expect_inherits(
  current = miloR_obj_approx,
  class = "miloR",
  info = "miloR class returned approx strategy"
)

expect_inherits(
  current = miloR_obj_bruteforce,
  class = "miloR",
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
  current = sum(get_differential_abundance_res(miloR_obj_index)$FDR <= 0.05) <
    3,
  info = "very few significant differences found in the correct case"
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

## exhaustive index ------------------------------------------------------------

# the exhaustive index scans every cell, so it returns the true nearest
# neighbour rather than an approximation
miloR_obj_exhaustive <- get_miloR_abundances_sc(
  object = sc_object,
  sample_id_col = "sample_id",
  miloR_params = params_sc_miloR(
    refinement_strategy = "index",
    index_type = "exhaustive"
  ),
  .verbose = FALSE
)

expect_inherits(
  current = miloR_obj_exhaustive,
  class = "miloR",
  info = "correct class returned for the exhaustive index"
)

expect_true(
  current = rs_set_similarity(
    as.character(get_index_cells(miloR_obj_exhaustive)),
    as.character(get_index_cells(miloR_obj_bruteforce)),
    TRUE
  ) >
    0.95,
  info = "the exhaustive index and brute force agree almost perfectly"
)

expect_error(
  current = params_sc_miloR(index_type = "nonsense"),
  info = "params_sc_miloR rejects an unknown index type"
)

## neighbourhood counting ------------------------------------------------------

# Rust returns the full grid, so unlike the old table() this cannot come back
# ragged
expect_equal(
  current = dim(miloR_obj_index$sample_counts),
  target = c(
    length(get_index_cells(miloR_obj_index)),
    nrow(unique(sc_object[["sample_id"]]))
  ),
  info = "the sample counts span every neighbourhood and every sample"
)

expect_true(
  current = all(miloR_obj_index$sample_counts >= 0),
  info = "the neighbourhood counts are non-negative"
)

# every neighbourhood holds its index cell plus its k neighbours, deduplicated,
# so its total across samples cannot exceed k + 1
expect_true(
  current = all(
    rowSums(miloR_obj_index$sample_counts) <= ncol(get_knn_mat(sc_object)) + 1
  ),
  info = "no cell is counted twice within a neighbourhood"
)

expect_equal(
  current = sum(miloR_obj_index$sample_counts),
  target = sum(miloR_obj_index$nhoods),
  info = "the counts total the non-zeros of the neighbourhood matrix"
)

## graph overlap ---------------------------------------------------------------

expect_equal(
  current = length(miloR_obj_index$nhood_overlap),
  target = length(get_index_cells(miloR_obj_index)),
  info = "one overlap value per neighbourhood"
)

expect_true(
  current = all(miloR_obj_index$nhood_overlap >= 0),
  info = "the overlaps are non-negative"
)

# the overlap is the row sums of t(nhoods) %*% nhoods with the diagonal zeroed,
# which the Rust never actually forms
intersect_mat <- Matrix::crossprod(miloR_obj_index$nhoods)
diag(intersect_mat) <- 0

expect_equal(
  current = miloR_obj_index$nhood_overlap,
  target = as.numeric(Matrix::rowSums(intersect_mat)),
  info = "the overlap matches the explicit crossproduct"
)

miloR_obj_overlap <- test_nhoods(
  x = miloR_obj_index_v2,
  design = ~grps,
  design_df = design_df,
  fdr_weighting = "graph-overlap"
)

expect_true(
  current = all(
    !is.na(get_differential_abundance_res(miloR_obj_overlap)$SpatialFDR)
  ),
  info = "graph-overlap weighting produces a spatial FDR"
)

expect_true(
  current = all(
    get_differential_abundance_res(miloR_obj_overlap)$SpatialFDR >= 0 &
      get_differential_abundance_res(miloR_obj_overlap)$SpatialFDR <= 1
  ),
  info = "the spatial FDR is a probability"
)

## spatial FDR -----------------------------------------------------------------

# a flat connectivity reduces the weighted step-up to plain Benjamini-Hochberg
set.seed(123L)
flat_p <- runif(200)

expect_equal(
  current = rs_spatial_fdr(
    p_values = flat_p,
    connectivity = rep(1, length(flat_p))
  ),
  target = stats::p.adjust(flat_p, method = "BH"),
  tolerance = 1e-12,
  info = "flat weights reduce the spatial FDR to Benjamini-Hochberg"
)

expect_true(
  current = all(is.na(rs_spatial_fdr(
    p_values = rep(NA_real_, 5),
    connectivity = rep(1, 5)
  ))),
  info = "non-finite p-values are carried through untouched"
)

expect_error(
  current = rs_spatial_fdr(p_values = flat_p, connectivity = c(1, 2)),
  info = "rs_spatial_fdr rejects mismatched lengths"
)

## coefficient selection -------------------------------------------------------

milo_by_name <- test_nhoods(
  x = miloR_obj_index_v2,
  design = ~grps,
  design_df = design_df,
  coef = "grpstrt"
)

expect_equal(
  current = get_differential_abundance_res(milo_by_name)$logFC,
  target = get_differential_abundance_res(miloR_obj_index_v2)$logFC,
  info = "naming the coefficient matches the default last-column choice"
)

milo_contrast <- test_nhoods(
  x = miloR_obj_index_v2,
  design = ~grps,
  design_df = design_df,
  contrast = c(0, 1)
)

expect_equal(
  current = get_differential_abundance_res(milo_contrast)$logFC,
  target = get_differential_abundance_res(miloR_obj_index_v2)$logFC,
  tolerance = 1e-8,
  info = "a unit contrast reproduces testing that coefficient directly"
)

## normalisation ---------------------------------------------------------------

milo_logms <- test_nhoods(
  x = miloR_obj_index_v2,
  design = ~grps,
  design_df = design_df,
  norm_method = "logMS"
)

expect_true(
  current = checkmate::testDataTable(
    get_differential_abundance_res(milo_logms)
  ),
  info = "logMS normalisation runs"
)

expect_error(
  current = test_nhoods(
    x = miloR_obj_index_v2,
    design = ~grps,
    design_df = design_df,
    norm_method = "nonsense"
  ),
  info = "test_nhoods rejects an unknown normalisation method"
)

## deprecated model getter -----------------------------------------------------

expect_warning(
  current = get_model_fit(miloR_obj_index),
  info = "get_model_fit is deprecated now that the fit lives in Rust"
)

expect_null(
  current = suppressWarnings(get_model_fit(miloR_obj_index)),
  info = "get_model_fit returns NULL"
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
  current = mean(miloR_obj_index$nhoods_info$majority_prop) > 0.9,
  info = "majority of the neighbourhoods are the same cell type"
)

## meld ------------------------------------------------------------------------

### distributions --------------------------------------------------------------

meld_sample_distributions <- sc_object[[c(
  "cell_id",
  "cell_grp",
  "sample_id_v2"
)]]

# sample 1, 2 are enriched for cell_type_1; 3, 4 for cell_type 2;
# sample 5, 6 for cell_type 3

unique_samples <- unique(meld_sample_distributions$sample_id_v2)

sample_to_cell <- lapply(unique_samples, FUN = function(x) {
  meld_sample_distributions[sample_id_v2 == x, cell_id]
})
names(sample_to_cell) <- unique_samples

### general tests --------------------------------------------------------------

meld_res <- meld_sc(
  object = sc_object,
  sample_id_col = "sample_id",
  .verbose = FALSE
)

expect_true(
  current = checkmate::testMatrix(
    meld_res$raw_scores,
    mode = "numeric",
    row.names = "named",
    col.names = "named"
  ),
  info = "meld results (raw) are a matrix"
)

expect_true(
  current = checkmate::testMatrix(
    meld_res$norm_scores,
    mode = "numeric",
    row.names = "named",
    col.names = "named"
  ),
  info = "meld results (norm) are a matrix"
)

expect_equal(
  current = dim(meld_res$raw_scores),
  target = c(nrow(meld_res$raw_scores), n_samples),
  info = "meld results has expected dimensionaliy"
)

expect_true(
  current = all(abs(rowSums(meld_res$norm_scores) - 1) < 1e-5),
  info = "meld normalised results are all close to 1"
)

### specific results -----------------------------------------------------------

meld_res_v2 <- meld_sc(
  object = sc_object,
  sample_id_col = "sample_id_v2",
  .verbose = FALSE
)

for (i in seq_along(sample_to_cell)) {
  sample_i <- names(sample_to_cell)[i]
  cells_i <- sample_to_cell[[i]]

  in_sample_i <- which(row.names(meld_res_v2$norm_scores) %in% cells_i)

  expect_true(
    current = mean(meld_res_v2$norm_scores[in_sample_i, sample_i]) >
      mean(meld_res_v2$norm_scores[-in_sample_i, sample_i]),
    info = sprintf("for %s the correct cells have higher meld scores", sample_i)
  )
}

# clean up ---------------------------------------------------------------------

sc_test_cleanup(test_temp_dir)
