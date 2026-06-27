# test data and params ---------------------------------------------------------

library(magrittr)

test_temp_dir <- file.path(
  tempdir(),
  "doublet_detection"
)

dir.create(test_temp_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot("Test directory does not exist" = dir.exists(test_temp_dir))

## parameters ------------------------------------------------------------------

n_doublets <- 200

## helper functions ------------------------------------------------------------

#' Helper function to calculate metrics
metrics_helper <- function(cm) {
  TP <- cm[2, 2]
  FP <- cm[1, 2]
  FN <- cm[2, 1]
  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  f1 <- 2 * (precision * recall) / (precision + recall)
  c(precision = precision, recall = recall, f1 = f1)
}

## initial data ----------------------------------------------------------------

syn_data <- generate_single_cell_test_data(seed = 123L)

ct1_idx <- which(syn_data$obs$cell_grp == "cell_type_1")
ct2_idx <- which(syn_data$obs$cell_grp == "cell_type_2")
ct3_idx <- which(syn_data$obs$cell_grp == "cell_type_3")

n_doublets_12 <- ceiling(n_doublets / 3)
n_doublets_23 <- ceiling(n_doublets / 3)
n_doublets_13 <- n_doublets - n_doublets_12 - n_doublets_23

doublet_pairs <- rbind(
  cbind(
    sample(ct1_idx, n_doublets_12, replace = TRUE),
    sample(ct2_idx, n_doublets_12, replace = TRUE)
  ),
  cbind(
    sample(ct2_idx, n_doublets_23, replace = TRUE),
    sample(ct3_idx, n_doublets_23, replace = TRUE)
  ),
  cbind(
    sample(ct1_idx, n_doublets_13, replace = TRUE),
    sample(ct3_idx, n_doublets_13, replace = TRUE)
  )
)

doublet_counts <- lapply(1:nrow(doublet_pairs), function(i) {
  syn_data$counts[doublet_pairs[i, 1], ] +
    syn_data$counts[doublet_pairs[i, 2], ]
})
doublet_matrix <- do.call(rbind, doublet_counts)
all_counts <- rbind(syn_data$counts, doublet_matrix)
n_total <- nrow(all_counts)

doublet_obs <- data.table::data.table(
  cell_id = sprintf("doublet_%04d", 1:n_doublets),
  cell_grp = "doublet",
  batch_index = 1,
  doublet = TRUE
)
new_obs <- data.table::rbindlist(list(
  syn_data$obs[, doublet := FALSE],
  doublet_obs
))
new_obs[, sample_id := rep(c("sample_A", "sample_B"), length.out = .N)]

## generate the object ---------------------------------------------------------

sc_object <- SingleCells(
  dir_data = test_temp_dir
)

# keep all cells for the sake of this
sc_object <- load_r_data(
  object = sc_object,
  counts = all_counts,
  obs = new_obs,
  var = syn_data$var,
  sc_qc_param = params_sc_min_quality(
    min_unique_genes = 0L,
    min_lib_size = 0L,
    min_cells = 0L
  ),
  streaming = 0L,
  .verbose = FALSE
)

# test scrublet ----------------------------------------------------------------

## rust logic ------------------------------------------------------------------

### optimal parameters according to the original paper -------------------------

optimal_params <- params_scrublet(
  normalisation = list(target_size = 1e4),
  pca = list(no_pcs = 15L),
  hvg = list(min_gene_var_pctl = 0.0),
  expected_doublet_rate = 0.2,
  sim_doublet_ratio = 1.0,
  n_bins = 100L
)

scrublet_res <- rs_sc_scrublet(
  f_path_gene = bixverse:::get_rust_count_gene_f_path(sc_object),
  f_path_cell = bixverse:::get_rust_count_cell_f_path(sc_object),
  cells_to_keep = get_cells_to_keep(sc_object),
  scrublet_params = optimal_params,
  seed = 42L,
  verbose = 0L,
  streaming = FALSE,
  return_combined_pca = TRUE,
  return_pairs = TRUE
)

expect_true(
  current = checkmate::qtest(
    scrublet_res$predicted_doublets,
    sprintf("B%s", nrow(new_obs))
  ),
  info = "rust scrublet: predicted doublets correct return type"
)

expect_true(
  current = checkmate::qtest(
    scrublet_res$doublet_scores_obs,
    sprintf("N%s", nrow(new_obs))
  ),
  info = "rust scrublet: doublet_scores_obs correct return type"
)

expect_true(
  current = checkmate::qtest(
    scrublet_res$doublet_scores_sim,
    sprintf("N%s", nrow(new_obs) * optimal_params$sim_doublet_ratio)
  ),
  info = "rust scrublet: doublet_scores_sim correct return type"
)

expect_true(
  current = checkmate::qtest(
    scrublet_res$doublet_errors_obs,
    sprintf("N%s", nrow(new_obs))
  ),
  info = "rust scrublet: doublet_scores_sim correct return type"
)

expect_true(
  current = length(scrublet_res$predicted_doublets) == nrow(all_counts),
  info = "rust scrublet: correct dimensions"
)

expect_true(
  current = checkmate::testMatrix(scrublet_res$pca, mode = "numeric"),
  info = "rust scrublet: matrix returned"
)

expect_true(
  current = checkmate::qtest(scrublet_res$pair_1, "I+"),
  info = "rust scrublet: parent 1 returned"
)

expect_true(
  current = checkmate::qtest(scrublet_res$pair_2, "I+"),
  info = "rust scrublet: parent 2 returned"
)

metrics <- metrics_helper(
  cm = table(
    new_obs$doublet,
    scrublet_res$predicted_doublets
  )
)

expect_true(
  current = metrics["recall"] >= 0.7,
  info = "rust scrublet: 'good' recall on synthetic data"
)

expect_true(
  current = metrics["f1"] >= 0.7,
  info = "rust scrublet: 'good' recall on synthetic data"
)

### other parameter versions ---------------------------------------------------

# log transform WITHOUT Z-scoring was seen to be optimal

params_full_norm <- params_scrublet(
  normalisation = list(
    target_size = 1e4,
    mean_center = TRUE,
    normalise_variance = TRUE
  ),
  pca = list(no_pcs = 15L),
  hvg = list(min_gene_var_pctl = 0.0),
  expected_doublet_rate = 0.2,
  sim_doublet_ratio = 1.0,
  n_bins = 100L
)

scrublet_res.full_norm <- rs_sc_scrublet(
  f_path_gene = bixverse:::get_rust_count_gene_f_path(sc_object),
  f_path_cell = bixverse:::get_rust_count_cell_f_path(sc_object),
  cells_to_keep = get_cells_to_keep(sc_object),
  scrublet_params = params_full_norm,
  seed = 42L,
  verbose = 0L,
  streaming = FALSE,
  return_combined_pca = FALSE,
  return_pairs = FALSE
)

metrics.full_norm <- metrics_helper(
  cm = table(
    new_obs$doublet,
    scrublet_res.full_norm$predicted_doublets
  )
)

# should still behave VERY similar

expect_true(
  current = metrics.full_norm["recall"] >= 0.7,
  info = "rust scrublet: 'good' recall on synthetic data"
)

expect_true(
  current = is.null(scrublet_res.full_norm$pca),
  info = "rust scrublet: PCA NOT returned"
)

expect_true(
  current = is.null(scrublet_res.full_norm$pair_1),
  info = "rust scrublet: parents NOT returned"
)

## s7 method -------------------------------------------------------------------

obj_res <- scrublet_sc(
  sc_object,
  scrublet_params = optimal_params,
  .verbose = FALSE
)

expect_true(
  current = checkmate::testClass(plot(obj_res), "ggplot"),
  info = "S7 scrublet: plot shown"
)

expect_true(
  current = checkmate::testClass(obj_res, "ScrubletRes"),
  info = "S7 scrublet: the correct class is being returned"
)

expect_true(
  current = checkmate::qtest(attr(obj_res, "cell_indices"), "I+"),
  info = "S7 scrublet: the cell indices are stored on the class"
)

expect_equivalent(
  current = obj_res$predicted_doublets,
  target = scrublet_res$predicted_doublets,
  info = "S7 scrublet: no weird changes during generation (called doublets)"
)

expect_true(
  current = cor(
    obj_res$doublet_scores_obs,
    scrublet_res$doublet_scores_obs,
    method = "spearman"
  ) >=
    0.99,
  info = "S7 scrublet: no weird changes during generation (obs scores)"
)

# update the thresholds
obj_res.up <- call_doublets_manual(obj_res, threshold = 0.175, .verbose = FALSE)

metrics_obj <- metrics_helper(
  cm = table(new_obs$doublet, obj_res.up$predicted_doublets)
)

expect_true(
  current = metrics_obj["recall"] >= 0.7,
  info = "rust scrublet: 'good' recall on synthetic data"
)

expect_true(
  current = metrics_obj["f1"] >= 0.7,
  info = "rust scrublet: 'good' recall on synthetic data"
)

obs_data <- get_data(obj_res.up)

expect_true(
  current = checkmate::testDataTable(obs_data),
  info = "getter on scrublet res working"
)

expect_true(
  current = checkmate::testNames(
    names(obs_data),
    must.include = c("doublet", "doublet_score", "cell_idx")
  ),
  info = "getter on scrublet res working - expected columns"
)

## s7 method with grouping ----------------------------------------------------

obj_res_grp <- scrublet_sc(
  sc_object,
  scrublet_params = optimal_params,
  group_by = "sample_id",
  .verbose = FALSE
)

expect_true(
  current = checkmate::testClass(obj_res_grp, "ScrubletRes"),
  info = "S7 scrublet grouped: correct class returned"
)

expect_true(
  current = isTRUE(attr(obj_res_grp, "grouped")),
  info = "S7 scrublet grouped: grouped attribute set"
)

expect_true(
  current = identical(attr(obj_res_grp, "group_by_col"), "sample_id"),
  info = "S7 scrublet grouped: group_by_col attribute set"
)

expect_true(
  current = checkmate::testNumeric(
    obj_res_grp$threshold,
    names = "named",
    len = 2L
  ),
  info = "S7 scrublet grouped: threshold is a named numeric vector"
)

expect_true(
  current = setequal(names(obj_res_grp$threshold), c("sample_A", "sample_B")),
  info = "S7 scrublet grouped: threshold has the expected sample names"
)

expect_true(
  current = checkmate::testCharacter(
    obj_res_grp$cell_groups,
    len = length(obj_res_grp$predicted_doublets)
  ),
  info = "S7 scrublet grouped: cell_groups aligns with cells"
)

expect_true(
  current = !is.unsorted(attr(obj_res_grp, "cell_indices")),
  info = "S7 scrublet grouped: cell_indices sorted ascending"
)

# call_doublets_manual on a specific sample
obj_res_grp.up <- call_doublets_manual(
  obj_res_grp,
  threshold = 0.2,
  for_sample = "sample_A",
  .verbose = FALSE
)

expect_true(
  current = obj_res_grp.up$threshold[["sample_A"]] == 0.2,
  info = "S7 scrublet grouped: threshold updated for the targeted sample"
)

expect_equivalent(
  current = obj_res_grp.up$threshold[["sample_B"]],
  target = obj_res_grp$threshold[["sample_B"]],
  info = "S7 scrublet grouped: other sample threshold untouched"
)

# default for_sample falls back to first
obj_res_grp.up2 <- call_doublets_manual(
  obj_res_grp,
  threshold = 0.2,
  .verbose = FALSE
)

expect_true(
  current = obj_res_grp.up2$threshold[[1L]] == 0.2,
  info = "S7 scrublet grouped: default for_sample updates first group"
)

expect_true(
  current = checkmate::testClass(plot(obj_res_grp), "ggplot"),
  info = "S7 scrublet grouped: plot works without for_sample"
)

expect_true(
  current = checkmate::testClass(
    plot(obj_res_grp, for_sample = "sample_B"),
    "ggplot"
  ),
  info = "S7 scrublet grouped: plot works with explicit for_sample"
)

obs_data_grp <- get_data(obj_res_grp)

expect_true(
  current = checkmate::testDataTable(obs_data_grp),
  info = "S7 scrublet grouped: get_data works"
)

expect_true(
  current = nrow(obs_data_grp) == length(obj_res_grp$predicted_doublets),
  info = "S7 scrublet grouped: get_data row count matches cells"
)

# test doublet detection -------------------------------------------------------

## rust logic ------------------------------------------------------------------

boost_params <- params_boost(
  hvg = list(min_gene_var_pctl = 0.0),
  pca = list(no_pcs = 10L),
  normalisation = list(target_size = 1e4),
  resolution = 0.5,
  voter_thresh = 0.25,
  n_iters = 10L
)

doublet_detection_res <- rs_sc_doublet_detection(
  f_path_gene = bixverse:::get_rust_count_gene_f_path(sc_object),
  f_path_cell = bixverse:::get_rust_count_cell_f_path(sc_object),
  cells_to_keep = get_cells_to_keep(sc_object),
  boost_params = boost_params,
  seed = 42L,
  verbose = 0L,
  streaming = FALSE
)

expect_true(
  current = checkmate::qtest(
    doublet_detection_res$doublet,
    sprintf("B%s", nrow(new_obs))
  ),
  info = paste(
    "rust boost classified doublet detection:",
    "predicted doublets correct return type"
  )
)

expect_true(
  current = checkmate::qtest(
    doublet_detection_res$doublet_score,
    sprintf("N%s", nrow(new_obs))
  ),
  info = paste(
    "rust boost classified doublet detection:",
    "doublet scores correct return type"
  )
)

expect_true(
  current = checkmate::qtest(
    doublet_detection_res$voting_avg,
    sprintf("N%s", nrow(new_obs))
  ),
  info = paste(
    "rust boost classified doublet detection:",
    "voting average correct return type"
  )
)

metrics <- metrics_helper(
  cm = table(
    new_obs$doublet,
    doublet_detection_res$doublet
  )
)

expect_true(
  current = metrics["recall"] >= 0.7,
  info = paste(
    "rust boost classified doublet detection:",
    "good recall"
  )
)

expect_true(
  current = metrics["f1"] >= 0.7,
  info = paste(
    "rust boost classified doublet detection:",
    "good f1 scores"
  )
)

## s7 method -------------------------------------------------------------------

obj_res <- doublet_detection_boost_sc(
  sc_object,
  boost_params = boost_params,
  .verbose = FALSE
)

expect_true(
  current = checkmate::testClass(obj_res, "BoostRes"),
  info = "S7 boost - correct class returned"
)

expect_equivalent(
  current = obj_res$doublet,
  target = doublet_detection_res$doublet,
  info = "S7 boost: no weird changes during generation (called doublets)"
)

obs_data <- get_data(obj_res)

expect_true(
  current = checkmate::testDataTable(obs_data),
  info = "getter on boost res working"
)

expect_true(
  current = checkmate::testNames(
    names(obs_data),
    must.include = c("doublet", "doublet_score", "cell_idx")
  ),
  info = "getter on boost res working - expected columns"
)

## s7 method with grouping ----------------------------------------------------

obj_res_grp <- doublet_detection_boost_sc(
  sc_object,
  boost_params = boost_params,
  group_by = "sample_id",
  .verbose = FALSE
)

expect_true(
  current = checkmate::testClass(obj_res_grp, "BoostRes"),
  info = "S7 boost grouped: correct class returned"
)

expect_true(
  current = isTRUE(attr(obj_res_grp, "grouped")),
  info = "S7 boost grouped: grouped attribute set"
)

expect_true(
  current = identical(attr(obj_res_grp, "group_by_col"), "sample_id"),
  info = "S7 boost grouped: group_by_col attribute set"
)

expect_true(
  current = checkmate::testCharacter(
    obj_res_grp$cell_groups,
    len = length(obj_res_grp$doublet)
  ),
  info = "S7 boost grouped: cell_groups aligns with cells"
)

expect_true(
  current = setequal(
    unique(obj_res_grp$cell_groups),
    c("sample_A", "sample_B")
  ),
  info = "S7 boost grouped: cell_groups contains expected samples"
)

expect_true(
  current = !is.unsorted(attr(obj_res_grp, "cell_indices")),
  info = "S7 boost grouped: cell_indices sorted ascending"
)

expect_true(
  current = length(obj_res_grp$doublet) == nrow(new_obs),
  info = "S7 boost grouped: all cells accounted for"
)

obs_data_grp <- get_data(obj_res_grp)

expect_true(
  current = checkmate::testDataTable(obs_data_grp),
  info = "S7 boost grouped: get_data works"
)

# test scdblfinder -------------------------------------------------------------

scdblfinder_params <- params_scdblfinder(
  pca = list(no_pcs = 10L),
  expected_doublet_rate = 0.17,
  n_genes = 50L,
  cxds_genes = 50L
)

## rust logic ------------------------------------------------------------------

scdblfinder_res <- rs_sc_scdblfinder(
  f_path_gene = bixverse:::get_rust_count_gene_f_path(sc_object),
  f_path_cell = bixverse:::get_rust_count_cell_f_path(sc_object),
  cell_indices = get_cells_to_keep(sc_object),
  params = scdblfinder_params,
  seed = 42L,
  return_features = FALSE,
  streaming = FALSE,
  verbose = 0L
)

expect_true(
  current = checkmate::qtest(
    scdblfinder_res$predicted_doublets,
    sprintf("B%s", nrow(new_obs))
  ),
  info = paste(
    "rust scdblfinder classified doublet detection:",
    "predicted doublets correct return type"
  )
)

expect_true(
  current = checkmate::qtest(
    scdblfinder_res$doublet_score,
    sprintf("N%s", nrow(new_obs))
  ),
  info = paste(
    "rust scdblfinder classified doublet detection:",
    "doublet scores correct return type"
  )
)

expect_true(
  current = checkmate::qtest(
    scdblfinder_res$threshold,
    "N1[0, 1]"
  ),
  info = paste(
    "rust scdblfinder classified doublet detection:",
    "threshold is correct type"
  )
)

metrics <- metrics_helper(
  cm = table(
    new_obs$doublet,
    scdblfinder_res$predicted_doublets
  )
)

expect_true(
  current = metrics["recall"] >= 0.7,
  info = paste(
    "rust scdblfinder classified doublet detection:",
    "good recall"
  )
)

expect_true(
  current = metrics["f1"] >= 0.7,
  info = paste(
    "rust scdblfinder classified doublet detection:",
    "good f1 scores"
  )
)

# check the other metrics

expect_true(
  current = mean(scdblfinder_res$weighted[1:1000]) <
    mean(scdblfinder_res$weighted[1001:1200]),
  info = "scdblfinder: the weighted scores are higher for doublets"
)

expect_true(
  current = mean(scdblfinder_res$cxds_scores[1:1000]) <
    mean(scdblfinder_res$cxds_scores[1001:1200]),
  info = "scdblfinder: the cxds scores are higher for doublets"
)

## s7 method -------------------------------------------------------------------

obj_res <- scdblfinder_sc(
  object = sc_object,
  scdblfinder_params = scdblfinder_params,
  .verbose = FALSE
)

expect_true(
  current = checkmate::testClass(obj_res, "ScDblFinderRes"),
  info = "S7 scDblFinder - correct class returned"
)

expect_equivalent(
  current = obj_res$predicted_doublets,
  target = scdblfinder_res$predicted_doublets,
  info = "S7 scDblFinder: no weird changes during generation (called doublets)"
)

obs_data <- get_data(obj_res)

expect_true(
  current = checkmate::testDataTable(obs_data),
  info = "getter on scdblfinder res working"
)

expect_true(
  current = checkmate::testNames(
    names(obs_data),
    must.include = c(
      "predicted_doublets",
      "doublet_score",
      "cxds_scores",
      "weighted",
      "cluster_labels",
      "cell_idx"
    )
  ),
  info = "getter on scdblfinder res working - expected columns"
)

expect_warning(
  current = get_feature_mat(obj_res),
  info = "get_features_mat() returns a warning"
)

obj_res <- scdblfinder_sc(
  object = sc_object,
  scdblfinder_params = scdblfinder_params,
  return_features = TRUE,
  .verbose = FALSE
)

expect_true(
  current = checkmate::testMatrix(
    get_feature_mat(obj_res),
    col.names = "named",
    row.names = "named",
    nrows = 1200L
  ),
  info = "get_features_mat() returns a feature matrix"
)

weighted_scores <- get_scores(obj_res, score_type = "weighted")

cxds_scores <- get_scores(obj_res, score_type = "cxds_scores")

expect_true(
  current = checkmate::qtest(weighted_scores, "N1200[0,1]"),
  info = "get_scores() returns a numeric of right length - weighted"
)

expect_true(
  current = checkmate::qtest(cxds_scores, "N1200[0,1]"),
  info = "get_scores() returns a numeric of right length - cxds"
)

## s7 method with grouping ----------------------------------------------------

obj_res_grp <- scdblfinder_sc(
  object = sc_object,
  scdblfinder_params = scdblfinder_params,
  group_by = "sample_id",
  .verbose = FALSE
)

expect_true(
  current = checkmate::testClass(obj_res_grp, "ScDblFinderRes"),
  info = "S7 scDblFinder grouped: correct class returned"
)

expect_true(
  current = isTRUE(attr(obj_res_grp, "grouped")),
  info = "S7 scDblFinder grouped: grouped attribute set"
)

expect_true(
  current = identical(attr(obj_res_grp, "group_by_col"), "sample_id"),
  info = "S7 scDblFinder grouped: group_by_col attribute set"
)

expect_true(
  current = checkmate::testNumeric(
    obj_res_grp$threshold,
    names = "named",
    len = 2L
  ),
  info = "S7 scDblFinder grouped: threshold is a named numeric vector"
)

expect_true(
  current = setequal(names(obj_res_grp$threshold), c("sample_A", "sample_B")),
  info = "S7 scDblFinder grouped: threshold has the expected sample names"
)

expect_true(
  current = checkmate::testCharacter(
    obj_res_grp$cell_groups,
    len = length(obj_res_grp$predicted_doublets)
  ),
  info = "S7 scDblFinder grouped: cell_groups aligns with cells"
)

expect_true(
  current = !is.unsorted(attr(obj_res_grp, "cell_indices")),
  info = "S7 scDblFinder grouped: cell_indices sorted ascending"
)

# cluster labels prefixed with sample name
expect_true(
  current = all(grepl("^(sample_A|sample_B)_", obj_res_grp$cluster_labels)),
  info = "S7 scDblFinder grouped: cluster labels prefixed with sample"
)

obs_data_grp <- get_data(obj_res_grp)

expect_true(
  current = checkmate::testDataTable(obs_data_grp),
  info = "S7 scDblFinder grouped: get_data works"
)

expect_true(
  current = checkmate::testNames(
    names(obs_data_grp),
    must.include = c(
      "predicted_doublets",
      "doublet_score",
      "cxds_scores",
      "weighted",
      "cluster_labels",
      "cell_idx"
    )
  ),
  info = "S7 scDblFinder grouped: get_data columns present"
)

# clean up ---------------------------------------------------------------------

on.exit(unlink(test_temp_dir, recursive = TRUE, force = TRUE), add = TRUE)
