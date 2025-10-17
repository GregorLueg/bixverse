# test data and params ---------------------------------------------------------

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

syn_data <- generate_single_cell_test_data()

# Get cell indices by cell type
ct1_idx <- which(syn_data$obs$cell_grp == "cell_type_1")
ct2_idx <- which(syn_data$obs$cell_grp == "cell_type_2")
ct3_idx <- which(syn_data$obs$cell_grp == "cell_type_3")

# Create cross-cell-type doublets (more detectable)
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

# update obs
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

## generate the object ---------------------------------------------------------

sc_object <- single_cell_exp(
  dir_data = tempdir()
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
  streaming = FALSE,
  .verbose = FALSE
)

# test scrublet ----------------------------------------------------------------

## rust logic ------------------------------------------------------------------

### optimal parameters according to the original paper -------------------------

# log transform WITHOUT Z-scoring was seen to be optimal
# see their page https://github.com/swolock/scrublet/blob/master/examples/demuxlet_example.ipynb

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
  verbose = FALSE,
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

# the data is so f--king weird it's better to
# use manual score here
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
  verbose = FALSE,
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
  current = metrics.full_norm["f1"] <= metrics["f1"],
  info = "rust scrublet: worse recall with bad paramters"
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
  current = checkmate::testClass(obj_res, "scrublet_res"),
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

expect_equivalent(
  current = obj_res$doublet_scores_obs,
  target = scrublet_res$doublet_scores_obs,
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

obs_data <- get_obs_data(obj_res.up)

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
  verbose = FALSE,
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
  current = checkmate::testClass(obj_res, "boost_res"),
  info = "S7 boost - correct class returned"
)

expect_equivalent(
  current = obj_res$doublet,
  target = doublet_detection_res$doublet,
  info = "S7 boost: no weird changes during generation (called doublets)"
)

obs_data <- get_obs_data(obj_res)

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
