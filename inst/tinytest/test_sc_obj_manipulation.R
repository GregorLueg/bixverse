# sc obj merging/slicing -------------------------------------------------------

source("helper_sc.R", local = TRUE)

library(magrittr)

test_temp_dir <- sc_test_dir("obj_manipulation")

## testing parameters ----------------------------------------------------------

# thresholds
min_lib_size <- 300L
min_genes_exp <- 45L
min_cells_exp <- 500L

## synthetic test data ---------------------------------------------------------

single_cell_test_data_1 <- sc_test_fixture(
  min_lib_size = min_lib_size,
  min_genes_exp = min_genes_exp,
  min_cells_exp = min_cells_exp,
  seed = 1L
)

single_cell_test_data_2 <- sc_test_fixture(
  min_lib_size = min_lib_size,
  min_genes_exp = min_genes_exp,
  min_cells_exp = min_cells_exp,
  seed = 2L
)

sc_qc_param <- sc_test_qc_params(single_cell_test_data_1, target_size = 1000)

path_obj_1 <- sc_test_dir("obj1", test_temp_dir)
path_obj_2 <- sc_test_dir("obj2", test_temp_dir)
path_merged <- sc_test_dir("merged", test_temp_dir)

## generate the two objects to merge -------------------------------------------

### object 1 -------------------------------------------------------------------

sc_object_1 <- sc_test_object(path_obj_1, single_cell_test_data_1)

features_1 <- get_sc_var(sc_object_1)[["gene_id"]]
cells_1 <- get_cell_names(sc_object_1)
counts_raw_1 <- sc_object_1[]
counts_norm_1 <- sc_object_1[,, assay = "norm"]

### object 2 -------------------------------------------------------------------

sc_object_2 <- sc_test_object(path_obj_2, single_cell_test_data_2)

features_2 <- get_sc_var(sc_object_2)[["gene_id"]]
cells_2 <- get_cell_names(sc_object_2)
counts_raw_2 <- sc_object_2[]
counts_norm_2 <- sc_object_2[,, assay = "norm"]

### merged object --------------------------------------------------------------

sc_object_merged <- SingleCells(dir_data = path_merged)

# tests ------------------------------------------------------------------------

## expected data ---------------------------------------------------------------

intersecting_features <- intersect(features_1, features_2)

## simple merge ----------------------------------------------------------------

sc_object_merged <- merge_sc_experiments(
  target = sc_object_merged,
  inputs = list(sc_object_1, sc_object_2),
  exp_ids = c("exp1", "exp2"),
  .verbose = FALSE
)

cells_merged <- unlist(sc_object_merged[["cell_id"]], use.names = FALSE)
features_merged <- get_sc_var(sc_object_merged)[["gene_id"]]

expect_equal(
  current = features_merged,
  target = intersecting_features,
  info = "merged object - right gene intersection"
)

expect_equal(
  current = cells_merged,
  target = c(sprintf("exp1_%s", cells_1), sprintf("exp2_%s", cells_2)),
  info = "merged object - right cell intersection"
)

expect_equivalent(
  current = sc_object_merged[],
  target = rbind(
    counts_raw_1[, intersecting_features],
    counts_raw_2[, intersecting_features]
  ),
  info = "merged object - raw counts correctly merged"
)

expect_equivalent(
  current = sc_object_merged[,, assay = "norm"],
  target = rbind(
    counts_norm_1[, intersecting_features],
    counts_norm_2[, intersecting_features]
  ),
  info = "merged object - norm counts correctly merged"
)

norm_counts_without_renormalisation <- sc_object_merged[,, assay = "norm"]

## merge with renormalisation --------------------------------------------------

sc_object_merged <- merge_sc_experiments(
  target = sc_object_merged,
  inputs = list(sc_object_1, sc_object_2),
  exp_ids = c("exp1", "exp2"),
  renormalise = TRUE,
  .verbose = FALSE
)

cells_merged <- unlist(sc_object_merged[["cell_id"]], use.names = FALSE)
features_merged <- get_sc_var(sc_object_merged)[["gene_id"]]

expect_equal(
  current = features_merged,
  target = intersecting_features,
  info = "merged object - right gene intersection"
)

expect_equal(
  current = cells_merged,
  target = c(sprintf("exp1_%s", cells_1), sprintf("exp2_%s", cells_2)),
  info = "merged object - right cell intersection"
)

expect_equivalent(
  current = sc_object_merged[],
  target = rbind(
    counts_raw_1[, intersecting_features],
    counts_raw_2[, intersecting_features]
  ),
  info = "merged object - raw counts correctly merged"
)

norm_counts_with_renormalisation <- sc_object_merged[,, assay = "norm"]

expect_true(
  current = any(
    norm_counts_with_renormalisation[1, ] !=
      norm_counts_without_renormalisation[1, ]
  ),
  info = "renormalisation should change some of the values (1)"
)

expect_true(
  current = any(
    norm_counts_with_renormalisation[100, ] !=
      norm_counts_without_renormalisation[100, ]
  ),
  info = "renormalisation should change some of the values (2)"
)

# clean up ---------------------------------------------------------------------

sc_test_cleanup(test_temp_dir)
