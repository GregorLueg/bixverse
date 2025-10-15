# sc processing ----------------------------------------------------------------

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

f_path_csr = file.path(tempdir(), "csr_test.h5ad")

obs_table <- single_cell_test_data$obs
obs_table$to_keep <- TRUE

write_h5ad_sc(
  f_path = f_path_csr,
  counts = single_cell_test_data$counts,
  obs = obs_table,
  var = single_cell_test_data$var,
  .verbose = FALSE
)

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

sc_qc_param = params_sc_min_quality(
  min_unique_genes = min_genes_exp,
  min_lib_size = min_lib_size,
  min_cells = min_cells_exp,
  target_size = 1000
)

## underlying class ------------------------------------------------------------

sc_object <- suppressWarnings(single_cell_exp(dir_data = tempdir()))

sc_object <- load_h5ad(
  object = sc_object,
  h5_path = path.expand(f_path_csr),
  sc_qc_param = sc_qc_param,
  .verbose = FALSE
)

# tests ------------------------------------------------------------------------

## function warnings -----------------------------------------------------------

expect_warning(
  current = find_hvg_sc(sc_object),
  info = "warning that no cells to keep were specified"
)

expect_warning(
  current = get_hvg(sc_object),
  info = "warning that hvgs can be found"
)

expect_warning(
  current = calculate_pca_sc(sc_object, no_pcs = 10L),
  info = "warning that no HVGs are detected"
)

expect_warning(
  current = find_neighbours_sc(sc_object),
  info = "warning that no PCA data are detected"
)

expect_warning(
  current = find_clusters_sc(sc_object),
  info = "warning that no kNN/sNN data was found"
)

## gene set proportions --------------------------------------------------------

gs_of_interest <- list(
  gs_1 = c("gene_001", "gene_002", "gene_003", "gene_004"),
  gs_2 = c("gene_096", "gene_097", "gene_100")
)

total_size <- Matrix::rowSums(counts_filtered)

props_gs_1 <- Matrix::rowSums(counts_filtered[, gs_of_interest$gs_1]) /
  total_size
props_gs_2 <- Matrix::rowSums(counts_filtered[, gs_of_interest$gs_2]) /
  total_size

sc_object <- gene_set_proportions_sc(
  sc_object,
  gs_of_interest,
  .verbose = FALSE
)

expect_equivalent(
  current = unlist(sc_object[["gs_1"]]),
  target = props_gs_1,
  tolerance = 10e-7,
  info = "gene proportion calculations gene set 1"
)

expect_equivalent(
  current = unlist(sc_object[["gs_2"]]),
  target = props_gs_2,
  tolerance = 10e-7,
  info = "gene proportion calculations gene set 2"
)

sc_object <- gene_set_proportions_sc(
  sc_object,
  gs_of_interest,
  .verbose = FALSE
)

expect_true(
  current = all(!c("gs_1.1", "gs_2.1") %in% colnames(sc_object[[]])),
  info = "overwriting of obs data works"
)

### streaming ------------------------------------------------------------------

sc_object <- gene_set_proportions_sc(
  sc_object,
  gs_of_interest,
  streaming = TRUE,
  .verbose = FALSE
)

expect_equivalent(
  current = unlist(sc_object[["gs_1"]]),
  target = props_gs_1,
  tolerance = 10e-7,
  info = "gene proportion calculations gene set 1 - streaming version"
)

expect_equivalent(
  current = unlist(sc_object[["gs_2"]]),
  target = props_gs_2,
  tolerance = 10e-7,
  info = "gene proportion calculations gene set 2 - streaming version"
)

## cells to keep logic ---------------------------------------------------------

threshold <- 0.05

cells_to_keep <- sc_object[[]][gs_2 < threshold, cell_id]

expect_true(
  current = length(cells_to_keep) > 600,
  info = "sensible cell filtering based on the threshold"
)

sc_object <- set_cell_to_keep(sc_object, cells_to_keep)

expect_true(
  current = all(
    unlist(sc_object[["cell_id"]]) == cells_to_keep
  ),
  info = "setting cells to keep removes them from the obs table"
)

cell_names_filtered <- get_cell_names(sc_object, filtered = TRUE)

expect_true(
  current = all(cell_names_filtered == cells_to_keep),
  info = "filter flag on cell names works"
)

counts_more_filtered <- counts_filtered[
  which(props_gs_2 < threshold),
]

expect_equivalent(
  current = get_cells_to_keep(sc_object),
  target = which(props_gs_2 < threshold) - 1, # zero indexed in Rust
  info = "cell to keep logic"
)

# the logic of retrieving cells will become more complicated here...

expect_equal(
  current = sc_object[,, use_cells_to_keep = TRUE],
  target = counts_more_filtered,
  info = "counts after setting cells to keep and using the parameter"
)

expect_equal(
  current = sc_object[,, use_cells_to_keep = FALSE],
  target = counts_filtered,
  info = "counts after setting cells to keep and NOT using the parameter"
)

## hvg selection ---------------------------------------------------------------

### vst version ----------------------------------------------------------------

#### r version -----------------------------------------------------------------

n_cells <- nrow(counts_more_filtered)

# calculate mean and variance
mean_values_r <- Matrix::colMeans(counts_more_filtered)
var_values_r <- matrixStats::colVars(as.matrix(counts_more_filtered)) *
  ((n_cells - 1) /
    n_cells)

valid_genes <- var_values_r > 0 & mean_values_r > 0

mean_log10 <- log10(mean_values_r[valid_genes])
var_log10 <- log10(var_values_r[valid_genes])

loess_fit <- loess(var_log10 ~ mean_log10, span = 0.3, degree = 2)
var_expected_log10 <- predict(loess_fit, mean_log10)
var_expected <- 10^var_expected_log10

clip_max <- sqrt(n_cells)

var_std <- sapply(1:ncol(counts_more_filtered), function(i) {
  gene_vals <- counts_more_filtered[, i]
  mean_i <- mean_values_r[i]
  sd_expected <- sqrt(var_expected[i])

  standardised <- (gene_vals - mean_i) / sd_expected

  standardised <- pmin(pmax(standardised, -clip_max), clip_max)

  # Variance of standardized values
  # population variance
  var(standardised) * ((n_cells - 1) / n_cells)
})

var_std[is.na(var_std)] <- 0

# need to 0-index
hvg_r <- order(var_std, decreasing = TRUE)[1:hvg_to_keep] - 1

#### rust part -----------------------------------------------------------------

##### direct load --------------------------------------------------------------

sc_object <- find_hvg_sc(
  object = sc_object,
  hvg_no = hvg_to_keep,
  .verbose = FALSE
)

var_data <- get_sc_var(sc_object)

expect_equivalent(
  current = var_data$mean,
  target = mean_values_r,
  tolerance = 10e-7,
  info = "Correct mean calculations for genes"
)

expect_equivalent(
  current = var_data$var,
  target = var_values_r,
  tolerance = 10e-7,
  info = "Correct variance calculations for genes"
)

expect_equivalent(
  current = var_data$var_std,
  target = var_std,
  tolerance = 10e-3,
  info = "Correct standardised variance calculations for genes"
)

hvg_rs <- get_hvg(sc_object)

expect_true(
  current = length(intersect(hvg_r, hvg_rs)) == hvg_to_keep,
  info = "Overlap in the detected HVGs"
)

##### streaming version --------------------------------------------------------

sc_object <- find_hvg_sc(
  object = sc_object,
  hvg_no = hvg_to_keep,
  .verbose = FALSE
)

var_data <- get_sc_var(sc_object)

expect_equivalent(
  current = var_data$mean,
  target = mean_values_r,
  tolerance = 10e-7,
  info = "Correct mean calculations for genes"
)

expect_equivalent(
  current = var_data$var,
  target = var_values_r,
  tolerance = 10e-7,
  info = "Correct variance calculations for genes"
)

expect_equivalent(
  current = var_data$var_std,
  target = var_std,
  tolerance = 10e-3,
  info = "Correct standardised variance calculations for genes"
)

hvg_rs <- get_hvg(sc_object)

expect_true(
  current = length(intersect(hvg_r, hvg_rs)) == hvg_to_keep,
  info = "Overlap in the detected HVGs"
)

## pca -------------------------------------------------------------------------

### r --------------------------------------------------------------------------

pca_input <- as.matrix(sc_object[,
  as.integer(get_hvg(sc_object) + 1),
  assay = "norm",
  return_format = "gene",
  use_cells_to_keep = TRUE
])

scaled_data <- scale(pca_input)

pca_r <- prcomp(pca_input, scale. = TRUE)

expected_names <- get_gene_names(sc_object)[hvg_r + 1]
actual_names <- colnames(pca_input)

### rust -----------------------------------------------------------------------

sc_object <- calculate_pca_sc(
  object = sc_object,
  no_pcs = no_pcs,
  randomised_svd = FALSE,
  .verbose = FALSE
)

expect_true(
  current = all.equal(
    abs(diag(cor(get_pca_factors(sc_object)[, 1:no_pcs], pca_r$x[, 1:no_pcs]))),
    rep(1, no_pcs),
    tolerance = 1e-8
  ),
  info = "PCA on single cell data compared to R"
)

sc_object <- calculate_pca_sc(
  object = sc_object,
  no_pcs = no_pcs,
  randomised_svd = TRUE,
  .verbose = FALSE
)

expect_true(
  current = all.equal(
    abs(diag(cor(get_pca_factors(sc_object)[, 1:no_pcs], pca_r$x[, 1:no_pcs]))),
    rep(1, no_pcs),
    tolerance = 1e-8
  ),
  info = "PCA on single cell data compared to R (randomised SVD)"
)

#### scaling within the rust function ------------------------------------------

# test if the scaling results in the same data
c(pca_factors, pca_loadings, pca_eigenvals, scaled) %<-%
  rs_sc_pca(
    f_path_gene = bixverse:::get_rust_count_gene_f_path(sc_object),
    no_pcs = no_pcs,
    random_svd = FALSE,
    cell_indices = get_cells_to_keep(sc_object),
    gene_indices = get_hvg(sc_object),
    seed = 42L,
    return_scaled = TRUE,
    verbose = FALSE
  )

expect_equivalent(current = scaled, target = scaled_data, tolerance = 1e-3)

## knn and snn -----------------------------------------------------------------

sc_object <- find_neighbours_sc(
  sc_object,
  neighbours_params = params_sc_neighbours(knn_algorithm = "hnsw"),
  .verbose = FALSE
)

expect_equal(
  current = dim(get_knn_mat(sc_object)),
  target = c(682, 15),
  info = "kNN matrix correctly returned"
)

expect_true(
  current = class(get_snn_graph(sc_object)) == "igraph",
  info = "igraph correctly returned"
)

## community detection ---------------------------------------------------------

sc_object <- find_clusters_sc(sc_object)


leiden_clusters <- unlist(sc_object[["leiden_clustering"]])
filter_ind <- is.na(leiden_clusters)

leiden_clusters <- leiden_clusters[!filter_ind]
cell_grps <- unlist(sc_object[["cell_grp"]])[!filter_ind]

f1_scores <- f1_score_confusion_mat(cell_grps, leiden_clusters)

expect_true(
  current = all(f1_scores > 0.95),
  info = "leiden clustering identifies the cell groups"
)

## dges ------------------------------------------------------------------------

### between two groups ---------------------------------------------------------

cell_names_1 <- sc_object[[]][leiden_clustering == 1, cell_id]
cell_names_2 <- sc_object[[]][leiden_clustering == 2, cell_id]

expect_error(
  current = find_markers_sc(
    object = sc_object,
    cells_1 = c(cell_names_1, "x"),
    cells_2 = cell_names_2
  ),
  info = "error if weird cells are selected"
)

dge_test <- find_markers_sc(
  object = sc_object,
  cells_1 = cell_names_1,
  cells_2 = cell_names_2,
  .verbose = FALSE
)

expected_upregulated <- sprintf("gene_%03d", 1:10)
expected_downregulated <- sprintf("gene_%03d", 21:30)

expect_true(
  current = all(
    expected_upregulated %in% dge_test[lfc > 0 & fdr <= 0.05, gene_id]
  ),
  info = "all expected up-regulated genes identified"
)

expect_true(
  current = all(
    expected_downregulated %in% dge_test[lfc < 0 & fdr <= 0.05, gene_id]
  ),
  info = "all expected down-regulated genes identified"
)

### find all markers -----------------------------------------------------------

dge_test_2 <- find_all_markers_sc(
  object = sc_object,
  column_of_interest = "leiden_clustering",
  .verbose = FALSE
)

all_cell_markers <- sprintf("gene_%03d", 1:30)

expect_true(
  current = all(
    all_cell_markers %in% dge_test_2[fdr <= 0.05, gene_id]
  ),
  info = "all expected cell markers identified"
)

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
  current = mean(auc_res_wilcox[
    "markers_cell_type_1",
    cells_per_cluster$cell_type_1
  ]) >=
    mean(auc_res_wilcox[
      "markers_cell_type_1",
      setdiff(
        get_cell_names(sc_object, filtered = TRUE),
        cells_per_cluster$cell_type_1
      )
    ]),
  info = paste(
    "auc values of expected cells",
    "with expected genes is higher (cell type 1)"
  )
)

expect_true(
  current = mean(auc_res_wilcox[
    "markers_cell_type_2",
    cells_per_cluster$cell_type_2
  ]) >=
    mean(auc_res_wilcox[
      "markers_cell_type_2",
      setdiff(
        get_cell_names(sc_object, filtered = TRUE),
        cells_per_cluster$cell_type_2
      )
    ]),
  info = paste(
    "auc values of expected cells",
    "with expected genes is higher (cell type 2)"
  )
)

expect_true(
  current = mean(auc_res_wilcox[
    "markers_cell_type_3",
    cells_per_cluster$cell_type_3
  ]) >=
    mean(auc_res_wilcox[
      "markers_cell_type_3",
      setdiff(
        get_cell_names(sc_object, filtered = TRUE),
        cells_per_cluster$cell_type_3
      )
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

## meta cell -------------------------------------------------------------------

meta_cell_data_v1 <- get_meta_cells_sc(
  object = sc_object,
  sc_meta_cell_params = params_sc_metacells(target_no_metacells = 10L),
  .verbose = FALSE
)

expect_equivalent(
  current = dim(meta_cell_data_v1$meta_raw_counts),
  target = c(10, 81),
  info = paste(
    "meta cell aggregation - correct dimensions raw counts"
  )
)

expect_equivalent(
  current = dim(meta_cell_data_v1$meta_norm_counts),
  target = c(10, 81),
  info = paste(
    "meta cell aggregation - correct dimensions norm counts"
  )
)

expect_true(
  current = checkmate::testClass(
    meta_cell_data_v1$meta_raw_counts,
    "dgRMatrix"
  ),
  info = paste("meta cell aggregation - expected return class")
)

expect_true(
  current = checkmate::testClass(
    meta_cell_data_v1$meta_norm_counts,
    "dgRMatrix"
  ),
  info = paste("meta cell aggregation - expected return class")
)

meta_cell_data_v2 <- get_meta_cells_sc(
  object = sc_object,
  sc_meta_cell_params = params_sc_metacells(target_no_metacells = 100L),
  .verbose = FALSE
)

expect_equivalent(
  current = dim(meta_cell_data_v2$meta_raw_counts),
  target = c(100, 81),
  info = paste(
    "meta cell aggregation - correct dimensions raw counts (second version)"
  )
)

expect_equivalent(
  current = dim(meta_cell_data_v2$meta_norm_counts),
  target = c(100, 81),
  info = paste(
    "meta cell aggregation - correct dimensions norm counts (second version)"
  )
)
