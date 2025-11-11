# sc processing ----------------------------------------------------------------

test_temp_dir <- file.path(
  tempdir(),
  paste0("test_", format(Sys.time(), "%Y%m%d_%H%M%S_"), sample(1000:9999, 1))
)
dir.create(test_temp_dir, recursive = TRUE)

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

counts_filtered <- single_cell_test_data$counts[cells_pass, genes_pass]

sc_qc_param <- params_sc_min_quality(
  min_unique_genes = min_genes_exp,
  min_lib_size = min_lib_size,
  min_cells = min_cells_exp,
  target_size = 1000
)

## underlying class ------------------------------------------------------------

sc_object <- single_cell_exp(dir_data = test_temp_dir)

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

object = sc_object
hotspot_params = params_sc_hotspot(
  model = "danb",
  knn = list(ann_dist = "cosine"),
  normalise = FALSE
)

rextendr::document()

hotspot_auto_cor <- rs_hotspot_autocor(
  f_path_genes = get_rust_count_gene_f_path(object),
  f_path_cells = get_rust_count_cell_f_path(object),
  embd = get_pca_factors(object),
  hotspot_params = hotspot_params,
  cells_to_keep = get_cells_to_keep(object),
  genes_to_use = get_gene_indices(
    object,
    gene_ids = get_gene_names(object),
    rust_index = TRUE
  ),
  streaming = FALSE,
  verbose = TRUE,
  seed = 42L
)
