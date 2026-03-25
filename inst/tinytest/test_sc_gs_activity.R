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
  current = suppressWarnings(aucell_sc(
    object = sc_object,
    gs_list = bad_list,
    auc_type = "auroc",
    .verbose = FALSE
  )),
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
  current = suppressWarnings(module_scores_sc(
    object = sc_object,
    gs_list = bad_list,
    .verbose = FALSE
  )),
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

available <- get_gene_names(sc_object)

vision_gs <- lapply(vision_gs, function(gs) {
  lapply(gs, intersect, available)
})

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

hotspot_gene_gene_cor <- generate_hotspot_membership(hotspot_gene_gene_cor)

hotspot_gene_gene_cor_streaming <- generate_hotspot_membership(
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

### synthetic data -------------------------------------------------------------

#### tfs -----------------------------------------------------------------------

# some of these are CLEARLY associated with other genes in terms of correlation
# structure
tfs <- sprintf("gene_%03i", c(1, 2, 11, 12, 21, 22, 50, 60, 70, 80, 90, 100))

tf_to_gene_map <- list(
  sprintf("gene_%03i", 1:10),
  sprintf("gene_%03i", 1:10),
  sprintf("gene_%03i", 11:20),
  sprintf("gene_%03i", 11:20),
  sprintf("gene_%03i", 21:30),
  sprintf("gene_%03i", 21:30)
) %>%
  `names<-`(sprintf("gene_%03i", c(1, 2, 11, 12, 21, 22)))

#### motif ranks ---------------------------------------------------------------

# Synthetic motif rankings: 100 genes x 6 motifs
# Lower rank = stronger regulatory potential
# motif_1: genes 1-5 ranked high (cell type A, first half)
# motif_2: genes 6-10 ranked high (cell type A, second half)
# motif_3: genes 11-15 ranked high (cell type B, first half)
# motif_4: genes 16-20 ranked high (cell type B, second half)
# motif_5: genes 21-25 ranked high (cell type C, first half)
# motif_6: genes 26-30 ranked high (cell type C, second half)

set.seed(42L)
n_genes <- 100L
n_motifs <- 6L
gene_names <- sprintf("gene_%03i", seq_len(n_genes))

# start with random rankings per motif
scenic_rankings <- matrix(0L, nrow = n_genes, ncol = n_motifs)
rownames(scenic_rankings) <- gene_names
colnames(scenic_rankings) <- sprintf("motif_%i", seq_len(n_motifs))

# each motif has 5 "true" target genes ranked in top 5, rest shuffled
motif_targets <- list(
  motif_1 = 1:5,
  motif_2 = 6:10,
  motif_3 = 11:15,
  motif_4 = 16:20,
  motif_5 = 21:25,
  motif_6 = 26:30
)

for (j in seq_len(n_motifs)) {
  targets <- motif_targets[[j]]
  non_targets <- setdiff(seq_len(n_genes), targets)
  # Targets get ranks 1-5, non-targets get shuffled 6-100
  scenic_rankings[targets, j] <- sample(seq_along(targets))
  scenic_rankings[non_targets, j] <- sample(
    (length(targets) + 1L):n_genes,
    length(non_targets)
  )
}

#### motif tf annotations ------------------------------------------------------

scenic_annot <- data.table::data.table(
  motif = c(
    "motif_1",
    "motif_2",
    "motif_3",
    "motif_4",
    "motif_5",
    "motif_6",
    "motif_1",
    "motif_3"
  ),
  TF = c(
    "gene_001",
    "gene_002",
    "gene_011",
    "gene_012",
    "gene_021",
    "gene_022",
    "gene_050",
    "gene_060"
  ),
  direct_annotation = c(
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    FALSE,
    FALSE
  ),
  inferred_orthology = c(
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    FALSE
  ),
  inferred_motif_sim = c(
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    TRUE,
    TRUE
  ),
  annotationSource = factor(c(
    rep("directAnnotation", 6),
    rep("inferredBy_MotifSimilarity", 2)
  )),
  description = "synthetic"
)

#### functions -----------------------------------------------------------------

check_tf_importances <- function(x, tf_to_gene_map = tf_to_gene_map) {
  # checks
  checkmate::assertMatrix(
    x,
    mode = "numeric",
    row.names = "named",
    col.names = "named"
  )
  checkmate::assertList(tf_to_gene_map)

  expected_imp <- rep(FALSE, times = length(tf_to_gene_map))
  represented_genes <- rownames(x)

  for (i in seq_along(tf_to_gene_map)) {
    tf_i <- names(tf_to_gene_map)[[i]]
    genes_i <- intersect(tf_to_gene_map[[i]], represented_genes)
    not_genes_i <- setdiff(represented_genes, genes_i)

    expected_imp[[i]] <- mean(x[genes_i, tf_i]) > mean(x[not_genes_i, tf_i])
  }

  names(expected_imp) <- names(tf_to_gene_map)

  expected_imp
}

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

importances_rf_batch_32 <- check_tf_importances(
  x = as.matrix(rf_scenic_res_batch_32),
  tf_to_gene_map = tf_to_gene_map
)

expect_true(
  current = all(importances_rf_batch_32),
  info = "scenic grn: rf (batch size 32) returns correct res"
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

importances_rf_batch_16 <- check_tf_importances(
  x = as.matrix(rf_scenic_res_batch_16),
  tf_to_gene_map = tf_to_gene_map
)

importances_rf_batch_1 <- check_tf_importances(
  x = as.matrix(rf_scenic_res_batch_1),
  tf_to_gene_map = tf_to_gene_map
)

expect_true(
  current = all(importances_rf_batch_16),
  info = "scenic grn: rf (batch size 16) returns correct res"
)

expect_true(
  current = all(importances_rf_batch_1),
  info = "scenic grn: rf (batch size 1) returns correct res"
)

##### add tf <> gene info ------------------------------------------------------

expect_warning(
  current = get_tf_to_gene(rf_scenic_res_batch_32),
  info = "scenic grn: no tf to gene added yet warning"
)

rf_scenic_res_batch_32 <- identify_tf_to_genes(
  rf_scenic_res_batch_32,
  k_tfs = 5L,
  k_genes = 10L,
  .verbose = FALSE
)

expect_true(
  current = checkmate::testDataTable(get_tf_to_gene(rf_scenic_res_batch_32)),
  info = "scenic grn: tf to gene added - right type"
)

expect_true(
  current = checkmate::testNames(
    names(get_tf_to_gene(rf_scenic_res_batch_32)),
    must.include = c("tf", "gene", "importance")
  ),
  info = "scenic grn: tf to gene added - columns"
)

nrow_before <- nrow(get_tf_to_gene(rf_scenic_res_batch_32))

rf_scenic_res_batch_32 <- tf_to_genes_correlations(
  x = rf_scenic_res_batch_32,
  object = sc_object,
  cor_filter = 0.05,
  .verbose = FALSE
)

expect_true(
  current = checkmate::testNames(
    names(get_tf_to_gene(rf_scenic_res_batch_32)),
    must.include = c("tf", "gene", "importance", "pairwise_cor")
  ),
  info = "scenic grn: tf to gene added - correlation added"
)

nrow_after <- nrow(get_tf_to_gene(rf_scenic_res_batch_32))

expect_true(
  current = nrow_after < nrow_before,
  info = "scenic grn: correlation filter behaves"
)

expect_true(
  current = all(
    get_tf_to_gene(rf_scenic_res_batch_32)[["pairwise_cor"]] >= 0.05
  ),
  info = "scenic grn: correlation filter behaves - test2"
)

rf_scenic_res_batch_32 <- tf_to_genes_motif_enrichment(
  x = rf_scenic_res_batch_32,
  motif_rankings = scenic_rankings,
  annot_data = scenic_annot,
  cis_target_params = params_cistarget(
    nes_threshold = 0.5,
    high_conf_cats = "directAnnotation",
    low_conf_cats = "inferredBy_MotifSimilarity"
  ),
  .verbose = FALSE
)

expect_true(
  current = checkmate::testDataTable(
    x = get_cistarget_res(rf_scenic_res_batch_32)
  ),
  info = "cistarget results added"
)

expect_true(
  current = checkmate::testNames(
    names(get_tf_to_gene(rf_scenic_res_batch_32)),
    must.include = c(
      "tf",
      "gene",
      "importance",
      "pairwise_cor",
      "in_leading_edge"
    )
  ),
  info = "scenic grn: motif - motif annotations added"
)

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

importances_et_batch_32 <- check_tf_importances(
  x = as.matrix(et_scenic_res_batch_32),
  tf_to_gene_map = tf_to_gene_map
)

importances_et_batch_16 <- check_tf_importances(
  x = as.matrix(et_scenic_res_batch_16),
  tf_to_gene_map = tf_to_gene_map
)

importances_et_batch_1 <- check_tf_importances(
  x = as.matrix(et_scenic_res_batch_1),
  tf_to_gene_map = tf_to_gene_map
)

expect_true(
  current = all(importances_et_batch_32),
  info = "scenic grn: rf (batch size 16) returns correct res"
)

expect_true(
  current = all(importances_et_batch_16),
  info = "scenic grn: rf (batch size 16) returns correct res"
)

expect_true(
  current = all(importances_et_batch_1),
  info = "scenic grn: rf (batch size 1) returns correct res"
)


#### boosted -------------------------------------------------------------------

grnboost2_scenic_res <- scenic_grn_sc(
  object = sc_object,
  tf_ids = tfs,
  genes_to_take = genes_to_keep_scenic,
  scenic_params = params_scenic(
    gene_batch_size = 1L,
    learner_type = "grnboost2"
  ),
  .verbose = FALSE
)

importances_grnboost <- check_tf_importances(
  x = as.matrix(grnboost2_scenic_res),
  tf_to_gene_map = tf_to_gene_map
)

expect_true(
  current = all(importances_grnboost),
  info = "scenic grn: grnboost returns correct res"
)

# clean up ---------------------------------------------------------------------

on.exit(unlink(test_temp_dir, recursive = TRUE, force = TRUE), add = TRUE)
