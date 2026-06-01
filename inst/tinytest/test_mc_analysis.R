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
  randomised_svd = TRUE,
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

## pca -------------------------------------------------------------------------

mc_object <- calculate_pca_sc(
  object = mc_object,
  no_pcs = no_pcs,
  randomised_svd = FALSE,
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
  randomised_svd = TRUE,
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

# clean up ---------------------------------------------------------------------

on.exit(unlink(test_temp_dir, recursive = TRUE, force = TRUE), add = TRUE)
