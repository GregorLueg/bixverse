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

write_h5ad_sc(
  f_path = f_path_csr,
  counts = single_cell_test_data$counts,
  obs = single_cell_test_data$obs,
  single_cell_test_data$var,
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
  current = all(unlist(sc_object[["cell_id"]]) == cells_to_keep),
  info = "setting genes to keep removes them from the obs table"
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
  as.integer(hvg_r + 1),
  assay = "norm",
  return_type = "gene",
  use_cells_to_keep = TRUE
])

pca_r <- prcomp(pca_input, scale. = TRUE)

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

cell_grps <- unlist(sc_object[["obs_cell_grp"]])
leiden_clusters <- unlist(sc_object[["leiden_clustering"]])

f1_scores <- f1_score_confusion_mat(cell_grps, leiden_clusters)

expect_true(
  current = all(f1_scores > 0.95),
  info = "leiden clustering identifies the cell groups"
)

cell_names_1 <- sc_object[[]][leiden_clustering == 1, cell_id]
cell_names_2 <- sc_object[[]][leiden_clustering == 2, cell_id]

?get_cell_indices

devtools::load_all()

rextendr::document()

get_cell_names(sc_object)

get_cell_indices(
  x = sc_object,
  cell_ids = cell_names_1,
  rust_index = TRUE
)


devtools::load_all()

x <- find_markers_sc(
  object = sc_object,
  cells_1 = c(cell_names_1),
  cells_2 = cell_names_2
)

get_rust_count_cell_f_path(sc_object)

rs_calculate_dge_mann_whitney(
  f_path = get_rust_count_cell_f_path(sc_object),
  cell_indices_1 = get_cell_indices(
    x = sc_object,
    cell_ids = cell_names_1,
    rust_index = TRUE
  ),
  cell_indices_2 = get_cell_indices(
    x = sc_object,
    cell_ids = cell_names_2,
    rust_index = TRUE
  ),
  min_prop = 0.6,
  alternative = "greater",
  verbose = TRUE
)
