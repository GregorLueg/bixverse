# meta cell processing and analysis --------------------------------------------

library(magrittr)

set.seed(123L)

test_temp_dir <- file.path(
  tempdir(),
  "mc_processing"
)

dir.create(test_temp_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot("Test directory does not exist" = dir.exists(test_temp_dir))

## testing parameters ----------------------------------------------------------

# thresholds
min_lib_size <- 300L
min_genes_exp <- 45L
min_cells_exp <- 500L
# meta cells
target_n_metacells <- 200L
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
  info = "mc processing - sensible amount of genes pass"
)

expect_true(
  current = length(cells_pass) > 800 & length(cells_pass) != 1000,
  info = "mc processing - sensible amount of cells pass"
)

## underlying single cell object -----------------------------------------------

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
  streaming = 0L,
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
  .verbose = FALSE
)

sc_object <- find_neighbours_sc(
  sc_object,
  neighbours_params = params_sc_neighbours(knn = list(k = 15L)),
  .verbose = FALSE
)

## generate hdWGCNA meta cells -------------------------------------------------

mc_object <- generate_bt_meta_cells_sc(
  sc_object,
  sc_meta_cell_params = params_sc_bt_metacells(
    target_no_metacells = target_n_metacells
  ),
  .verbose = FALSE
)

# dominant original cell type per meta cell - used by AUCell assertions
original_cell_type <- unlist(sc_object[["cell_grp"]])
mc_dominant_type <- purrr::map_chr(
  mc_object[[]]$original_cell_idx,
  function(idx) {
    types <- original_cell_type[idx]
    names(which.max(table(types)))
  }
)

# tests ------------------------------------------------------------------------

## hvg -------------------------------------------------------------------------

expect_warning(
  current = calculate_pca_sc(mc_object, no_pcs = no_pcs, .verbose = FALSE),
  info = "mc - warning when no HVGs identified yet"
)

### vst ------------------------------------------------------------------------

mc_object <- find_hvg_sc(
  object = mc_object,
  hvg_no = hvg_to_keep,
  .verbose = FALSE
)

mc_hvg <- get_hvg(mc_object)

expect_true(
  current = checkmate::qtest(mc_hvg, sprintf("I%i", hvg_to_keep)),
  info = "mc hvg - correct length of indices returned"
)

# var table is updated with hvg metrics
mc_var <- get_sc_var(mc_object)

expect_true(
  current = all(c("mean", "var", "var_std") %in% names(mc_var)),
  info = "mc hvg - var table populated with hvg metrics"
)

# the signal genes (1-30) should dominate the HVGs - it's the same signal
# carried into the meta cells
expect_true(
  current = length(intersect(mc_hvg, 1:30)) >=
    length(intersect(mc_hvg, 31:100)),
  info = "mc hvg - bulk of signal genes recovered as HVGs"
)

### dispersion -----------------------------------------------------------------

mc_object <- find_hvg_sc(
  object = mc_object,
  hvg_no = hvg_to_keep,
  hvg_params = params_sc_hvg(method = "dispersion"),
  .verbose = FALSE
)

mc_var_disp <- get_sc_var(mc_object)

expect_true(
  current = all(
    c("mean", "dispersion", "dispersion_scaled", "bin") %in% names(mc_var_disp)
  ),
  info = "mc hvg dispersion - var table populated with dispersion metrics"
)

mc_hvg_disp <- get_hvg(mc_object)

expect_true(
  current = checkmate::qtest(mc_hvg_disp, sprintf("I%i", hvg_to_keep)),
  info = "mc hvg dispersion - correct number of HVGs"
)

expect_true(
  current = length(intersect(mc_hvg_disp, 1:30)) >=
    length(intersect(mc_hvg_disp, 31:100)),
  info = "mc hvg dispersion - signal genes dominate the HVGs"
)

### meanvarbin -----------------------------------------------------------------

mc_object <- find_hvg_sc(
  object = mc_object,
  hvg_no = hvg_to_keep,
  hvg_params = params_sc_hvg(method = "meanvarbin"),
  .verbose = FALSE
)

mc_var_mvb <- get_sc_var(mc_object)

expect_true(
  current = all(
    c("mean", "dispersion", "dispersion_scaled", "bin") %in% names(mc_var_mvb)
  ),
  info = "mc hvg meanvarbin - var table populated with metrics"
)

mc_hvg_mvb <- get_hvg(mc_object)

expect_true(
  current = checkmate::qtest(mc_hvg_mvb, sprintf("I%i", hvg_to_keep)),
  info = "mc hvg meanvarbin - correct number of HVGs"
)

expect_true(
  current = length(intersect(mc_hvg_mvb, 1:30)) >=
    length(intersect(mc_hvg_mvb, 31:100)),
  info = "mc hvg meanvarbin - signal genes dominate the HVGs"
)

### get_hvg_data_sc ------------------------------------------------------------

# reset to vst for downstream tests
mc_object <- find_hvg_sc(
  object = mc_object,
  hvg_no = hvg_to_keep,
  .verbose = FALSE
)

mc_current_hvg <- get_hvg(mc_object)
mc_current_var <- get_sc_var(mc_object)

#### all meta cells ------------------------------------------------------------

mc_hvg_dt_full <- get_hvg_data_sc(
  object = mc_object,
  hvg_no = hvg_to_keep,
  .verbose = FALSE
)

expect_true(
  current = checkmate::testDataTable(mc_hvg_dt_full),
  info = "mc get_hvg_data_sc - returns a data.table"
)

expect_true(
  current = all(
    c(
      "gene_idx",
      "gene_id",
      "mean",
      "var",
      "var_std",
      "is_hvg",
      "hvg_rank"
    ) %in%
      names(mc_hvg_dt_full)
  ),
  info = "mc get_hvg_data_sc - expected columns present"
)

expect_true(
  current = sum(mc_hvg_dt_full$is_hvg) == hvg_to_keep,
  info = "mc get_hvg_data_sc - correct number of HVGs flagged"
)

expect_true(
  current = all(!is.na(mc_hvg_dt_full[is_hvg == TRUE, hvg_rank])) &
    all(is.na(mc_hvg_dt_full[is_hvg == FALSE, hvg_rank])),
  info = "mc get_hvg_data_sc - hvg_rank NA exactly when is_hvg is FALSE"
)

# cell_ids = NULL should mirror the state-mutating version on all meta cells
expect_true(
  current = length(intersect(
    mc_hvg_dt_full[is_hvg == TRUE, gene_idx],
    mc_current_hvg
  )) ==
    hvg_to_keep,
  info = "mc get_hvg_data_sc - NULL cell_ids matches find_hvg_sc"
)

#### biological subset ---------------------------------------------------------

ct1_mc_ids <- mc_object[[]]$meta_cell_id[mc_dominant_type == "cell_type_1"]

expect_true(
  current = length(ct1_mc_ids) > 10,
  info = "mc get_hvg_data_sc - sensible number of cell_type_1 meta cells"
)

mc_hvg_dt_subset <- get_hvg_data_sc(
  object = mc_object,
  cell_ids = ct1_mc_ids,
  hvg_no = hvg_to_keep,
  .verbose = FALSE
)

expect_true(
  current = sum(mc_hvg_dt_subset$is_hvg) == hvg_to_keep,
  info = "mc get_hvg_data_sc - subset returns correct number of HVGs"
)

expect_false(
  current = isTRUE(all.equal(mc_hvg_dt_subset$mean, mc_hvg_dt_full$mean)),
  info = "mc get_hvg_data_sc - subset stats differ from full stats"
)

#### state preservation --------------------------------------------------------

expect_equal(
  current = get_hvg(mc_object),
  target = mc_current_hvg,
  info = "mc get_hvg_data_sc - object HVGs unchanged after call"
)

expect_equivalent(
  current = get_sc_var(mc_object)$var_std,
  target = mc_current_var$var_std,
  info = "mc get_hvg_data_sc - object var table unchanged after call"
)

## pca -------------------------------------------------------------------------

mc_object <- calculate_pca_sc(
  object = mc_object,
  no_pcs = no_pcs,
  pca_params = params_sc_pca(randomised = FALSE),
  .verbose = FALSE
)

pca_factors <- get_pca_factors(mc_object)
pca_loadings <- get_pca_loadings(mc_object)
pca_singular <- get_pca_singular_val(mc_object)

expect_true(
  current = checkmate::testMatrix(
    pca_factors,
    mode = "numeric",
    nrows = target_n_metacells,
    ncols = no_pcs,
    row.names = "named",
    col.names = "named"
  ),
  info = "mc pca - factors correct shape and naming"
)

expect_true(
  current = checkmate::testMatrix(
    pca_loadings,
    mode = "numeric",
    nrows = hvg_to_keep,
    ncols = no_pcs,
    row.names = "named",
    col.names = "named"
  ),
  info = "mc pca - loadings correct shape and naming"
)

expect_true(
  current = checkmate::qtest(pca_singular, sprintf("N%i", no_pcs)),
  info = "mc pca - singular values correct length"
)

# randomised SVD path
mc_object <- calculate_pca_sc(
  object = mc_object,
  no_pcs = no_pcs,
  pca_params = params_sc_pca(),
  .verbose = FALSE
)

expect_true(
  current = checkmate::testMatrix(
    get_pca_factors(mc_object),
    mode = "numeric",
    nrows = target_n_metacells,
    ncols = no_pcs
  ),
  info = "mc pca - randomised SVD path returns correct shape"
)

pca_factors_randomised <- get_pca_factors(mc_object)

expect_true(
  current = all(abs(diag(cor(pca_factors, pca_factors_randomised))) >= 0.99),
  info = "randomised svd and svd return similar values"
)

## knn / snn -------------------------------------------------------------------

knn_k <- 5L

mc_object <- find_neighbours_sc(
  mc_object,
  neighbours_params = params_sc_neighbours(knn = list(k = knn_k)),
  .verbose = FALSE
)

expect_equal(
  current = dim(get_knn_mat(mc_object)),
  target = c(target_n_metacells, knn_k),
  info = "mc knn - correct shape of kNN matrix"
)

expect_equal(
  current = dim(get_knn_dist(mc_object)),
  target = c(target_n_metacells, knn_k),
  info = "mc knn - correct shape of kNN distance matrix"
)

expect_true(
  current = class(get_snn_graph(mc_object)) == "igraph",
  info = "mc snn - igraph returned"
)

expect_equal(
  current = igraph::vcount(get_snn_graph(mc_object)),
  target = target_n_metacells,
  info = "mc snn - graph vertex count matches meta cell count"
)

## aucell ----------------------------------------------------------------------

auc_gene_sets <- list(
  markers_cell_type_1 = sprintf("gene_%03d", 1:10),
  markers_cell_type_2 = sprintf("gene_%03d", 11:20),
  markers_cell_type_3 = sprintf("gene_%03d", 21:30)
)

bad_list <- list(markers = sample(letters, 10))

expect_error(
  current = suppressWarnings(aucell_sc(
    object = mc_object,
    gs_list = bad_list,
    auc_type = "auroc",
    .verbose = FALSE
  )),
  info = "mc aucell - error when no gene set entries match"
)

auc_res_wilcox <- aucell_sc(
  object = mc_object,
  gs_list = auc_gene_sets,
  auc_type = "wilcox",
  .verbose = FALSE
)

auc_res_auroc <- aucell_sc(
  object = mc_object,
  gs_list = auc_gene_sets,
  auc_type = "auroc",
  .verbose = FALSE
)

expect_true(
  current = checkmate::testMatrix(
    auc_res_wilcox,
    mode = "numeric",
    nrows = target_n_metacells,
    ncols = length(auc_gene_sets),
    row.names = "named",
    col.names = "named"
  ),
  info = "mc aucell wilcox - correct matrix shape and naming"
)

expect_true(
  current = checkmate::testMatrix(
    auc_res_auroc,
    mode = "numeric",
    nrows = target_n_metacells,
    ncols = length(auc_gene_sets),
    row.names = "named",
    col.names = "named"
  ),
  info = "mc aucell auroc - correct matrix shape and naming"
)

# meta cells dominated by a given cell type should score higher on the
# matching marker set
mc_ids <- mc_object[[]]$meta_cell_id

for (ct in c("cell_type_1", "cell_type_2", "cell_type_3")) {
  mc_in <- mc_ids[mc_dominant_type == ct]
  mc_out <- mc_ids[mc_dominant_type != ct]
  marker_col <- sprintf("markers_%s", ct)

  expect_true(
    current = mean(auc_res_wilcox[mc_in, marker_col]) >
      mean(auc_res_wilcox[mc_out, marker_col]),
    info = sprintf(
      "mc aucell wilcox - %s markers higher in matching meta cells",
      ct
    )
  )

  expect_true(
    current = mean(auc_res_auroc[mc_in, marker_col]) >
      mean(auc_res_auroc[mc_out, marker_col]),
    info = sprintf(
      "mc aucell auroc - %s markers higher in matching meta cells",
      ct
    )
  )
}

expect_true(
  current = all(diag(cor(auc_res_wilcox, auc_res_auroc)) >= 0.99),
  info = "mc aucell - wilcox and auroc highly correlated"
)

## scenic ----------------------------------------------------------------------

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

check_tf_importances <- function(x, tf_to_gene_map) {
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

### gene filter ----------------------------------------------------------------

genes_to_keep_scenic <- scenic_gene_filter_sc(
  object = mc_object,
  scenic_params = params_scenic(min_counts = 500L),
  .verbose = FALSE
)

expect_true(
  current = checkmate::qtest(genes_to_keep_scenic, "S+"),
  info = "mc scenic - gene filter returns character vector"
)

expect_true(
  current = length(intersect(tfs, genes_to_keep_scenic)) >= 6,
  info = "mc scenic - sensible number of TFs survive gene filter"
)

### tf <> gene importances -----------------------------------------------------

#### random forest -------------------------------------------------------------

rf_scenic_res <- scenic_grn_sc(
  object = mc_object,
  tf_ids = tfs,
  genes_to_take = genes_to_keep_scenic,
  scenic_params = params_scenic(gene_batch_size = 32L),
  .verbose = FALSE
)

expect_true(
  current = checkmate::testClass(rf_scenic_res, "ScenicGrn"),
  info = "mc scenic rf - correct class returned"
)

expect_true(
  current = checkmate::testMatrix(
    as.matrix(rf_scenic_res),
    mode = "numeric",
    row.names = "named",
    col.names = "named"
  ),
  info = "mc scenic rf - as.matrix() returns correctly named numeric matrix"
)

importances_rf <- check_tf_importances(
  x = as.matrix(rf_scenic_res),
  tf_to_gene_map = tf_to_gene_map
)

expect_true(
  current = all(importances_rf),
  info = "mc scenic rf - expected TF -> gene importances recovered"
)

#### extra trees ---------------------------------------------------------------

et_scenic_res <- scenic_grn_sc(
  object = mc_object,
  tf_ids = tfs,
  genes_to_take = genes_to_keep_scenic,
  scenic_params = params_scenic(
    gene_batch_size = 32L,
    learner_type = "extratrees"
  ),
  .verbose = FALSE
)

importances_et <- check_tf_importances(
  x = as.matrix(et_scenic_res),
  tf_to_gene_map = tf_to_gene_map
)

expect_true(
  current = all(importances_et),
  info = "mc scenic extratrees - expected TF -> gene importances recovered"
)

#### grnboost2 -----------------------------------------------------------------

grnboost_scenic_res <- scenic_grn_sc(
  object = mc_object,
  tf_ids = tfs,
  genes_to_take = genes_to_keep_scenic,
  scenic_params = params_scenic(
    gene_batch_size = 32L,
    learner_type = "grnboost2"
  ),
  .verbose = FALSE
)

importances_grnboost <- check_tf_importances(
  x = as.matrix(grnboost_scenic_res),
  tf_to_gene_map = tf_to_gene_map
)

expect_true(
  current = all(importances_grnboost),
  info = "mc scenic grnboost2 - expected TF -> gene importances recovered"
)

## nmf -------------------------------------------------------------------------

### helpers --------------------------------------------------------------------

# For each NMF component, find which true cell type has the highest mean H
# activation. Returns a vector of length k mapping component index to cell type.
.nmf_component_to_celltype <- function(h, cell_grp_per_col) {
  apply(h, 1, function(activations) {
    means <- tapply(activations, cell_grp_per_col, mean)
    names(means)[which.max(means)]
  })
}

# For each NMF component, find which signal gene block (1:10, 11:20, 21:30)
# has the highest mean W loading. Returns a vector mapping component to block.
.nmf_component_to_gene_block <- function(w) {
  gene_blocks <- list(
    cell_type_1 = sprintf("gene_%03d", 1:10),
    cell_type_2 = sprintf("gene_%03d", 11:20),
    cell_type_3 = sprintf("gene_%03d", 21:30)
  )
  available <- rownames(w)
  apply(w, 2, function(loadings) {
    block_means <- sapply(gene_blocks, function(genes) {
      mean(loadings[intersect(genes, available)])
    })
    names(which.max(block_means))
  })
}

### tests ----------------------------------------------------------------------

mc_nmf_k <- 3L

### nmf_sc on MetaCells --------------------------------------------------------

mc_nmf_res <- nmf_sc(
  object = mc_object,
  k = mc_nmf_k,
  preprocessing = "none",
  use_second_layer = TRUE,
  .verbose = FALSE
)

expect_true(
  current = checkmate::testClass(mc_nmf_res, "NmfResult"),
  info = "mc nmf_sc: NmfResult class returned"
)

mc_hvg_n <- length(get_hvg(mc_object))

expect_equal(
  current = dim(get_w(mc_nmf_res)),
  target = c(mc_hvg_n, mc_nmf_k),
  info = "mc nmf_sc: W has shape (n_hvg, k)"
)

expect_equal(
  current = dim(get_h(mc_nmf_res)),
  target = c(mc_nmf_k, target_n_metacells),
  info = "mc nmf_sc: H has shape (k, n_metacells)"
)

expect_true(
  current = all(get_w(mc_nmf_res) >= 0) & all(get_h(mc_nmf_res) >= 0),
  info = "mc nmf_sc: W and H are non-negative"
)

#### signal recovery via dominant cell type -----------------------------------

mc_h <- get_h(mc_nmf_res)
mc_w <- get_w(mc_nmf_res)

mc_grp_per_col <- mc_dominant_type[match(colnames(mc_h), mc_ids)]

mc_comp_to_ct <- .nmf_component_to_celltype(mc_h, mc_grp_per_col)
mc_comp_to_block <- .nmf_component_to_gene_block(mc_w)

expect_equal(
  current = sort(unique(mc_comp_to_ct)),
  target = c("cell_type_1", "cell_type_2", "cell_type_3"),
  info = "mc nmf_sc: H components bijectively cover the three cell types"
)

expect_equal(
  current = sort(unique(mc_comp_to_block)),
  target = c("cell_type_1", "cell_type_2", "cell_type_3"),
  info = "mc nmf_sc: W components bijectively cover the three gene blocks"
)

expect_equal(
  current = mc_comp_to_ct,
  target = mc_comp_to_block,
  info = "mc nmf_sc: per-component, H cell type matches W gene block"
)

#### get_data slots into obs ---------------------------------------------------

mc_nmf_obs_dt <- get_data(mc_nmf_res)

expect_true(
  current = checkmate::testDataTable(mc_nmf_obs_dt, nrows = target_n_metacells),
  info = "mc nmf_sc: get_data returns data.table of correct shape"
)

expect_true(
  current = isTRUE(attr(mc_nmf_obs_dt, "is_obs")),
  info = "mc nmf_sc: get_data carries the is_obs attribute"
)

expect_true(
  current = all(sprintf("comp_%02d", 1:mc_nmf_k) %in% names(mc_nmf_obs_dt)),
  info = "mc nmf_sc: get_data has one column per component"
)

### stabilised_nmf_sc on MetaCells ---------------------------------------------

mc_n_runs <- 4L

mc_stab_res <- stabilised_nmf_sc(
  object = mc_object,
  k = mc_nmf_k,
  preprocessing = "none",
  use_second_layer = TRUE,
  n_runs = mc_n_runs,
  .verbose = FALSE
)

expect_true(
  current = checkmate::testClass(mc_stab_res, "StabilisedNmfResult"),
  info = "mc stabilised_nmf_sc: StabilisedNmfResult class returned"
)

expect_equal(
  current = dim(get_w(mc_stab_res)),
  target = c(mc_hvg_n, mc_nmf_k * mc_n_runs),
  info = "mc stabilised_nmf_sc: w_all has shape (n_hvg, k * n_runs)"
)

expect_true(
  current = is.list(get_h(mc_stab_res)) &
    length(get_h(mc_stab_res)) == mc_n_runs,
  info = "mc stabilised_nmf_sc: h_per_run is a list of length n_runs"
)

expect_equal(
  current = mc_stab_res$best_idx,
  target = which.min(mc_stab_res$losses),
  info = "mc stabilised_nmf_sc: best_idx matches minimum loss run"
)

#### get_best_run round-trip ---------------------------------------------------

mc_best <- get_best_run(mc_stab_res)

expect_true(
  current = checkmate::testClass(mc_best, "NmfResult"),
  info = "mc get_best_run: returns an NmfResult"
)

mc_k_idx <- ((mc_stab_res$best_idx - 1L) *
  mc_nmf_k +
  1L):(mc_stab_res$best_idx *
  mc_nmf_k)

expect_equivalent(
  current = unname(get_w(mc_stab_res)[, mc_k_idx]),
  target = unname(get_w(mc_best)),
  info = "mc get_best_run: W matches the corresponding slice of w_all"
)

expect_equivalent(
  current = unname(get_h(mc_stab_res)[[mc_stab_res$best_idx]]),
  target = unname(get_h(mc_best)),
  info = "mc get_best_run: H matches the corresponding entry in h_per_run"
)

mc_best_w <- get_w(mc_best)
mc_best_h <- get_h(mc_best)

mc_best_grp_per_col <- mc_dominant_type[match(colnames(mc_best_h), mc_ids)]
mc_best_comp_to_ct <- .nmf_component_to_celltype(mc_best_h, mc_best_grp_per_col)
mc_best_comp_to_block <- .nmf_component_to_gene_block(mc_best_w)

expect_equal(
  current = sort(unique(mc_best_comp_to_ct)),
  target = c("cell_type_1", "cell_type_2", "cell_type_3"),
  info = "mc stabilised_nmf_sc: best run H covers all cell types"
)

expect_equal(
  current = mc_best_comp_to_ct,
  target = mc_best_comp_to_block,
  info = "mc stabilised_nmf_sc: best run components consistent across W and H"
)

# clean up ---------------------------------------------------------------------

on.exit(unlink(test_temp_dir, recursive = TRUE, force = TRUE), add = TRUE)
