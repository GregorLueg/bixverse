# sc processing ----------------------------------------------------------------

library(magrittr)

test_temp_dir <- file.path(
  tempdir(),
  "gs_activitiy"
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

## underlying class ------------------------------------------------------------

sc_object <- SingleCells(dir_data = test_temp_dir)

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

sc_object <- calculate_pca_sc(sc_object, no_pcs = no_pcs, .verbose = FALSE)

# tests ------------------------------------------------------------------------

## aucell ----------------------------------------------------------------------

auc_gene_sets <- list(
  markers_cell_type_1 = sprintf("gene_%03d", 1:10),
  markers_cell_type_2 = sprintf("gene_%03d", 11:20),
  markers_cell_type_3 = sprintf("gene_%03d", 21:30)
)

bad_list <- list(markers = sample(letters, 10))

auc_res_wilcox <- aucell_sc(
  object = sc_object,
  gs_list = auc_gene_sets,
  auc_type = "wilcox",
  .verbose = FALSE
)

auc_res_auroc <- aucell_sc(
  object = sc_object,
  gs_list = auc_gene_sets,
  auc_type = "auroc",
  .verbose = FALSE
)

obs_table_red <- sc_object[[c("cell_id", "cell_grp")]]

cells_per_cluster <- split(
  obs_table_red$cell_id,
  obs_table_red$cell_grp
)

expect_error(
  current = aucell_sc(
    object = sc_object,
    gs_list = bad_list,
    auc_type = "auroc",
    .verbose = FALSE
  ),
  info = paste("aucell: error when provided a list where nothing matches.")
)

expect_true(
  current = checkmate::testMatrix(
    auc_res_wilcox,
    mode = "numeric",
    ncols = length(auc_gene_sets),
    nrows = length(get_cells_to_keep(sc_object)),
    col.names = "named",
    row.names = "named"
  ),
  info = paste(
    "aucell results wilcox statistic what you'd expect"
  )
)

expect_true(
  current = checkmate::testMatrix(
    auc_res_auroc,
    mode = "numeric",
    ncols = length(auc_gene_sets),
    nrows = length(get_cells_to_keep(sc_object)),
    col.names = "named",
    row.names = "named"
  ),
  info = paste(
    "aucell results auroc statistic what you'd expect"
  )
)

expect_true(
  current = mean(auc_res_wilcox[
    cells_per_cluster$cell_type_1,
    "markers_cell_type_1"
  ]) >=
    mean(auc_res_wilcox[
      setdiff(
        get_cell_names(sc_object, filtered = TRUE),
        cells_per_cluster$cell_type_1
      ),
      "markers_cell_type_1"
    ]),
  info = paste(
    "auc values of expected cells",
    "with expected genes is higher (cell type 1)"
  )
)

expect_true(
  current = mean(auc_res_wilcox[
    cells_per_cluster$cell_type_2,
    "markers_cell_type_2"
  ]) >=
    mean(auc_res_wilcox[
      setdiff(
        get_cell_names(sc_object, filtered = TRUE),
        cells_per_cluster$cell_type_2
      ),
      "markers_cell_type_2"
    ]),
  info = paste(
    "auc values of expected cells",
    "with expected genes is higher (cell type 2)"
  )
)

expect_true(
  current = mean(auc_res_wilcox[
    cells_per_cluster$cell_type_3,
    "markers_cell_type_3"
  ]) >=
    mean(auc_res_wilcox[
      setdiff(
        get_cell_names(sc_object, filtered = TRUE),
        cells_per_cluster$cell_type_3
      ),
      "markers_cell_type_3"
    ]),
  info = paste(
    "auc values of expected cells",
    "with expected genes is higher (cell type 3)"
  )
)

expect_true(
  current = all(diag(cor(auc_res_wilcox, auc_res_auroc)) >= 0.99),
  info = paste(
    "auc values between the two methods are highly correlated"
  )
)


## module scores ---------------------------------------------------------------

module_scores <- module_scores_sc(
  object = sc_object,
  gs_list = auc_gene_sets,
  .verbose = FALSE
)

module_scores_streaming <- module_scores_sc(
  object = sc_object,
  gs_list = auc_gene_sets,
  .verbose = FALSE,
  streaming = TRUE
)

expect_error(
  current = module_scores_sc(
    object = sc_object,
    gs_list = bad_list,
    .verbose = FALSE
  ),
  info = paste(
    "module-scores: error when provided a list where nothing matches."
  )
)

expect_true(
  current = checkmate::testMatrix(
    module_scores,
    mode = "numeric",
    ncols = length(auc_gene_sets),
    nrows = length(get_cells_to_keep(sc_object)),
    col.names = "named",
    row.names = "named"
  ),
  info = paste(
    "module scores what you'd expect"
  )
)

expect_true(
  current = checkmate::testMatrix(
    module_scores_streaming,
    mode = "numeric",
    ncols = length(auc_gene_sets),
    nrows = length(get_cells_to_keep(sc_object)),
    col.names = "named",
    row.names = "named"
  ),
  info = paste(
    "module scores what you'd expect (streaming version)"
  )
)

expect_equivalent(
  current = module_scores,
  target = module_scores_streaming,
  info = paste(
    "the two versions are equivalent"
  )
)

expect_true(
  current = mean(module_scores[
    cells_per_cluster$cell_type_1,
    "markers_cell_type_1"
  ]) >=
    mean(module_scores[
      setdiff(
        get_cell_names(sc_object, filtered = TRUE),
        cells_per_cluster$cell_type_1
      ),
      "markers_cell_type_1"
    ]),
  info = paste(
    "modules score values of expected cells",
    "with expected genes is higher (cell type 1)"
  )
)

expect_true(
  current = mean(module_scores[
    cells_per_cluster$cell_type_2,
    "markers_cell_type_2"
  ]) >=
    mean(module_scores[
      setdiff(
        get_cell_names(sc_object, filtered = TRUE),
        cells_per_cluster$cell_type_2
      ),
      "markers_cell_type_2"
    ]),
  info = paste(
    "modules score values of expected cells",
    "with expected genes is higher (cell type 2)"
  )
)

expect_true(
  current = mean(module_scores[
    cells_per_cluster$cell_type_3,
    "markers_cell_type_3"
  ]) >=
    mean(module_scores[
      setdiff(
        get_cell_names(sc_object, filtered = TRUE),
        cells_per_cluster$cell_type_3
      ),
      "markers_cell_type_3"
    ]),
  info = paste(
    "modules score values of expected cells",
    "with expected genes is higher (cell type 3)"
  )
)

## vision pathway score --------------------------------------------------------

set.seed(42L)

vision_gs <- list(
  gs_1 = list(
    pos = sprintf("gene_%03d", 1:10),
    neg = sprintf("gene_%03d", 11:30)
  ),
  gs_2 = list(
    pos = sprintf("gene_%03d", 11:20),
    neg = sprintf("gene_%03d", c(1:10, 21:40))
  ),
  gs_3 = list(
    pos = sprintf("gene_%03d", 21:30)
  ),
  gs_4 = list(
    pos = sample(get_gene_names(sc_object)[30:80], 11),
    neg = sample(get_gene_names(sc_object)[30:80], 9)
  ),
  gs_5 = list(
    pos = sample(get_gene_names(sc_object)[30:80], 6),
    neg = sample(get_gene_names(sc_object)[30:80], 14)
  ),
  gs_6 = list(
    pos = sample(get_gene_names(sc_object)[30:80], 2),
    neg = sample(get_gene_names(sc_object)[30:80], 8)
  )
)

vision_res <- vision_sc(
  object = sc_object,
  gs_list = vision_gs,
  .verbose = FALSE
)

expect_true(
  current = checkmate::testMatrix(
    vision_res,
    mode = "numeric",
    ncols = length(vision_gs),
    nrows = length(get_cells_to_keep(sc_object)),
    col.names = "named",
    row.names = "named"
  ),
  info = paste(
    "vision (matrix) outputs what you'd expect"
  )
)

expect_true(
  current = sign(mean(vision_res[, 1])) == -1,
  info = "vision - delta gene set v1 has negative avg scores"
)

expect_true(
  current = sign(mean(vision_res[, 2])) == -1,
  info = "vision - delta gene set v2 has negative avg scores"
)

expect_true(
  current = sign(mean(vision_res[, 3])) == 1,
  info = "vision - gene set withou delta has positive signs"
)

### vision with auto-correlation -----------------------------------------------

vision_res_auto <- vision_w_autocor_sc(
  object = sc_object,
  gs_list = vision_gs,
  embd_to_use = "pca",
  vision_params = params_sc_vision(n_perm = 100L),
  .verbose = FALSE
)

expect_equivalent(
  current = vision_res_auto$vision_matrix,
  target = vision_res,
  info = "vision w autocor - matrix is still sensible"
)

expect_true(
  current = checkmate::testDataTable(
    vision_res_auto$auto_cor_dt,
    nrows = length(vision_gs)
  ),
  info = "vision w autocor - autocor data.table is correct"
)

expect_equivalent(
  current = vision_res_auto$auto_cor_dt$p_val < 0.05,
  target = c(rep(TRUE, 3), rep(FALSE, 3)),
  info = "vision w autocor - correct gene sets are significant"
)

vision_res_auto_v2 <- vision_w_autocor_sc(
  object = sc_object,
  gs_list = vision_gs[1:4],
  embd_to_use = "pca",
  vision_params = params_sc_vision(n_perm = 100L),
  .verbose = FALSE
)

expect_true(
  current = checkmate::testMatrix(
    vision_res_auto_v2$vision_matrix,
    mode = "numeric",
    ncols = 4, # reduced
    nrows = length(get_cells_to_keep(sc_object)),
    col.names = "named",
    row.names = "named"
  ),
  info = paste(
    "vision (matrix) outputs with 4 gene sets what you'd expect"
  )
)

## hotspots --------------------------------------------------------------------

### auto-correlations ----------------------------------------------------------

hotspot_autocor_danb_res <- hotspot_autocor_sc(
  object = sc_object,
  .verbose = FALSE
)

expect_true(
  current = checkmate::testDataTable(hotspot_autocor_danb_res, nrows = 81L),
  info = "hotspot - local correlation expected outputs"
)

expect_true(
  current = all(hotspot_autocor_danb_res$fdr[1:30] <= 0.05),
  info = "hotspot - local correlations significant were expected - DANB"
)

expect_true(
  current = mean(hotspot_autocor_danb_res$gaerys_c[1:30]) >=
    mean(hotspot_autocor_danb_res$gaerys_c[31:nrow(hotspot_autocor_danb_res)]),
  info = "hotspot - average Gaery's C higher in expected genes - DANB"
)

hotspot_autocor_normal_res <- hotspot_autocor_sc(
  object = sc_object,
  hotspot_params = params_sc_hotspot(model = "normal"),
  .verbose = FALSE
)

expect_true(
  current = all(hotspot_autocor_normal_res$fdr[1:30] <= 0.05),
  info = "hotspot - local correlations significant were expected - normal"
)

expect_true(
  current = mean(hotspot_autocor_normal_res$gaerys_c[1:30]) >=
    mean(hotspot_autocor_normal_res$gaerys_c[
      31:nrow(hotspot_autocor_normal_res)
    ]),
  info = "hotspot - average Gaery's C higher in expected genes - normal"
)

#### streaming -----------------------------------------------------------------

hotspot_autocor_danb_res_streaming <- hotspot_autocor_sc(
  object = sc_object,
  streaming = TRUE,
  .verbose = FALSE
)

expect_equal(
  current = hotspot_autocor_danb_res_streaming,
  target = hotspot_autocor_danb_res,
  info = "hotspot - streaming engine behaves for auto-correlations"
)

#### gene subset ---------------------------------------------------------------

hotspot_autocor_danb_res_small <- hotspot_autocor_sc(
  object = sc_object,
  genes_to_take = get_gene_names(sc_object)[1:50],
  .verbose = FALSE
)

expect_true(
  current = checkmate::testDataTable(
    hotspot_autocor_danb_res_small,
    nrows = 50L
  ),
  info = "hotspot - local correlation expected outputs - reduced genes"
)

#### cell subset ---------------------------------------------------------------

cells_to_take <- sc_object[[]][cell_grp != "cell_type_1", cell_id]

hotspot_autocor_danb_res_cell_subset <- hotspot_autocor_sc(
  object = sc_object,
  cells_to_take = cells_to_take,
  .verbose = FALSE
)

# the first 10 genes should have much lower spatial correlation and Z-scores

expect_true(
  current = mean(hotspot_autocor_danb_res$gaerys_c[1:10]) >
    mean(hotspot_autocor_danb_res_cell_subset$gaerys_c[1:10]),
  info = "hotspot - local corr: sample removal as expected impact on gaery's c"
)

expect_true(
  current = mean(hotspot_autocor_danb_res$z_score[1:10]) >
    mean(hotspot_autocor_danb_res_cell_subset$z_score[1:10]),
  info = "hotspot - local corr: sample removal as expected impact on z scores"
)

### gene-gene correlations -----------------------------------------------------

hotspot_gene_gene_cor <- hotspot_gene_cor_sc(
  object = sc_object,
  genes_to_take = hotspot_autocor_danb_res[fdr <= 0.05, gene_id],
  .verbose = FALSE
)

hotspot_gene_gene_cor_streaming <- hotspot_gene_cor_sc(
  object = sc_object,
  genes_to_take = hotspot_autocor_danb_res[fdr <= 0.05, gene_id],
  streaming = TRUE,
  .verbose = FALSE
)

expect_equal(
  current = hotspot_gene_gene_cor$z,
  target = hotspot_gene_gene_cor_streaming$z,
  # tiny differences due to the batch-wise calculations here
  tolerance = 1e-5,
  info = paste(
    "hotspot - gene-gene corr:",
    "streaming and in-memory z matrix equivalence"
  )
)

expect_equal(
  current = hotspot_gene_gene_cor$cor,
  target = hotspot_gene_gene_cor_streaming$cor,
  # tiny differences due to the batch-wise calculations here
  tolerance = 1e-5,
  info = paste(
    "hotspot - gene-gene corr:",
    "streaming and in-memory cor matrix equivalence"
  )
)

hotspot_gene_gene_cor <- set_hotspot_membership(hotspot_gene_gene_cor)

hotspot_gene_gene_cor_streaming <- set_hotspot_membership(
  hotspot_gene_gene_cor_streaming
)

expect_equivalent(
  current = get_hotspot_membership(hotspot_gene_gene_cor)$cluster_member,
  target = get_hotspot_membership(
    hotspot_gene_gene_cor_streaming
  )$cluster_member,
  info = paste(
    "hotspot - gene-gene corr:",
    "streaming and in-memory cluster membership equivalence"
  )
)

expect_true(
  current = sum(
    !is.na(hotspot_gene_gene_cor$module_memership$cluster_member)
  ) ==
    30L,
  info = paste(
    "hotspot - gene-gene corr:",
    "expected cluster membership"
  )
)

## scenic ----------------------------------------------------------------------

### gene filtering -------------------------------------------------------------

genes_pre_filter <- get_sc_var(sc_object)

filtered_genes <- genes_pre_filter[no_cells_exp >= 500L, gene_id]

genes_to_keep_scenic <- scenic_gene_filter_sc(
  object = sc_object,
  # absurdly high, but releted to synthetic data
  scenic_params = params_scenic(min_counts = 500L),
  .verbose = FALSE
)

expect_true(
  current = checkmate::qtest(genes_to_keep_scenic, "S+"),
  info = "scenic_gene_filter returning right type"
)

expect_true(
  current = length(setdiff(filtered_genes, genes_to_keep_scenic)) == 0,
  info = "scenic_gene_filter returning the right genes"
)

### scenic tf <> gene ----------------------------------------------------------

tfs <- sprintf("gene_%03i", c(1, 2, 11, 12, 21, 22, 50, 60, 70, 80, 90, 100))

expect_true(
  current = length(intersect(tfs, filtered_genes)) > 6,
  info = "sensible tf numbers pass"
)

#### rf ------------------------------------------------------------------------

# test various batch sizes and check that stuff behaves
rf_scenic_res_batch_32 <- scenic_grn_sc(
  object = sc_object,
  tf_ids = tfs,
  genes_to_take = genes_to_keep_scenic,
  scenic_params = params_scenic(gene_batch_size = 32L),
  .verbose = FALSE
)


str(rf_scenic_res_batch_32)

expect_true(
  current = checkmate::testClass(rf_scenic_res_batch_32, "ScenicGrn"),
  info = "scenic correct class returned"
)

expect_true(
  current = checkmate::testMatrix(
    as.matrix(rf_scenic_res_batch_32),
    mode = "numeric",
    row.names = "named",
    col.names = "named"
  ),
  info = "scenic correct as.matrix() behaves"
)

##### other batch sizes --------------------------------------------------------

rf_scenic_res_batch_16 <- scenic_grn_sc(
  object = sc_object,
  tf_ids = tfs,
  genes_to_take = genes_to_keep_scenic,
  scenic_params = params_scenic(gene_batch_size = 16L),
  .verbose = FALSE
)

rf_scenic_res_batch_1 <- scenic_grn_sc(
  object = sc_object,
  tf_ids = tfs,
  genes_to_take = genes_to_keep_scenic,
  scenic_params = params_scenic(gene_batch_size = 1L),
  .verbose = FALSE
)

##### add tf <> gene info ------------------------------------------------------

#### extra trees ---------------------------------------------------------------

et_scenic_res_batch_32 <- scenic_grn_sc(
  object = sc_object,
  tf_ids = tfs,
  genes_to_take = genes_to_keep_scenic,
  scenic_params = params_scenic(
    gene_batch_size = 32L,
    learner_type = "extratrees"
  ),
  .verbose = FALSE
)

et_scenic_res_batch_16 <- scenic_grn_sc(
  object = sc_object,
  tf_ids = tfs,
  genes_to_take = genes_to_keep_scenic,
  scenic_params = params_scenic(
    gene_batch_size = 16L,
    learner_type = "extratrees"
  ),
  .verbose = FALSE
)

et_scenic_res_batch_1 <- scenic_grn_sc(
  object = sc_object,
  tf_ids = tfs,
  genes_to_take = genes_to_keep_scenic,
  scenic_params = params_scenic(
    gene_batch_size = 1L,
    learner_type = "extratrees"
  ),
  .verbose = FALSE
)

#### boosted -------------------------------------------------------------------

grnboost2_scenic_res_batch_1 <- scenic_grn_sc(
  object = sc_object,
  tf_ids = tfs,
  genes_to_take = genes_to_keep_scenic,
  scenic_params = params_scenic(
    gene_batch_size = 1L,
    learner_type = "grnboost2"
  ),
  .verbose = FALSE
)

# clean up ---------------------------------------------------------------------

on.exit(unlink(test_temp_dir, recursive = TRUE, force = TRUE), add = TRUE)


# dasda

class(x)

x = rf_scenic_res_batch_32
k_tfs = 5L
k_genes = 10L
min_importance = NULL

gene_tf_imp <- x$importance_matrix

# check if users wants to go from gene -> TF
tf_to_gene_1 <- if (!is.null(k_tfs)) {
  gene_to_tf <- rs_top_k_targets(
    matrix = gene_tf_imp,
    k = k_tfs,
    margin = 1L,
    min_value = min_importance
  )

  tapply(
    rep(names(gene_to_tf), lengths(gene_to_tf)),
    unlist(gene_to_tf),
    unique
  )
} else {
  NULL
}

# check if users wants to go from TF -> gene
tf_to_gene_2 <- if (!is.null(k_tfs)) {
  rs_top_k_targets(
    matrix = gene_tf_imp,
    k = k_tfs,
    margin = 2L,
    min_value = min_importance
  )
} else {
  NULL
}


# combine lists
tf_to_gene_final <- if (!is.null(tf_to_gene_1) & !is.null(tf_to_gene_2)) {
  purrr::map2(tf_to_gene_1, tf_to_gene_2, \(tg_1, tg_2) {
    sort(union(tg_1, tg_2))
  })
} else if (!is.null(tf_to_gene_1)) {
  tf_to_gene_1
} else {
  tf_to_gene_2
}

tf_gene_dt <- data.table::data.table(
  tf = rep(names(tf_to_gene_final), lengths(tf_to_gene_final)),
  target_gene = unlist(tf_to_gene_final, use.names = FALSE)
)


get_gene_indices(x = sc_object, gene_ids = tfs, rust_index = TRUE)
