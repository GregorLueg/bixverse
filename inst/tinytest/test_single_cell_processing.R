# sc processing ----------------------------------------------------------------

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

# thresholds
# absurd numbers, but this is due to the synthetic data
min_lib_size <- 300L
min_genes_exp <- 45L
min_cells_exp <- 500L

genes_pass <- which(
  Matrix::colSums(single_cell_test_data$counts != 0) >= min_cells_exp
)

cells_pass <- which(
  (Matrix::rowSums(single_cell_test_data$counts[, genes_pass]) >=
    min_lib_size) &
    (Matrix::rowSums(single_cell_test_data$counts[, genes_pass] != 0) >=
      min_genes_exp)
)

# test 1
expect_true(
  current = length(genes_pass) > 80 & length(genes_pass) != 100,
  info = "sc processing - sensible amount of genes pass"
)

# test 2
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

# test 3
expect_warning(
  current = find_hvg_sc(sc_object),
  info = "warning that no cells to keep were specified"
)

# test 4
expect_warning(
  current = calculate_pca_sc(sc_object, no_pcs = 10L),
  info = "warning that no HVGs are detected"
)

# test 5
expect_warning(
  current = find_neighbours_sc(sc_object),
  info = "warning that no PCA data are detected"
)

# test 6
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

# test 7
expect_equivalent(
  current = unlist(sc_object[["gs_1"]]),
  target = props_gs_1,
  tolerance = 10e-7,
  info = "gene proportion calculations gene set 1"
)

# test 8
expect_equivalent(
  current = unlist(sc_object[["gs_2"]]),
  target = props_gs_2,
  tolerance = 10e-7,
  info = "gene proportion calculations gene set 2"
)

# rerun and test that columns are not duplicated
sc_object <- gene_set_proportions_sc(
  sc_object,
  gs_of_interest,
  .verbose = FALSE
)

# test 9
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

# test 10
expect_equivalent(
  current = unlist(sc_object[["gs_1"]]),
  target = props_gs_1,
  tolerance = 10e-7,
  info = "gene proportion calculations gene set 1 - streaming version"
)

# test 11
expect_equivalent(
  current = unlist(sc_object[["gs_2"]]),
  target = props_gs_2,
  tolerance = 10e-7,
  info = "gene proportion calculations gene set 2 - streaming version"
)

## cells to keep logic ---------------------------------------------------------

threshold <- 0.05

cells_to_keep <- sc_object[[]][gs_2 < threshold, cell_id]

# test 12
expect_true(
  current = length(cells_to_keep) > 600,
  info = "sensible cell filtering based on the threshold"
)

sc_object <- set_cell_to_keep(sc_object, cells_to_keep)

# test 13
expect_true(
  current = all(unlist(sc_object[["cell_id"]]) == cells_to_keep),
  info = "setting genes to keep removes them from the obs table"
)

cell_names_filtered <- get_cell_names(sc_object, filtered = TRUE)

# test 14
expect_true(
  current = all(cell_names_filtered == cells_to_keep),
  info = "filter flag on cell names works"
)

counts_more_filtered <- counts_filtered[
  which(props_gs_2 < threshold),
]

# test 15
expect_equivalent(
  current = get_cells_to_keep(sc_object),
  target = which(props_gs_2 < threshold) - 1, # zero indexed in Rust
  info = "cell to keep logic"
)

# the logic of retrieving cells will become more complicated here...

# test 16
expect_equal(
  current = sc_object[,, use_cells_to_keep = TRUE],
  target = counts_more_filtered,
  info = "counts after setting cells to keep and using the parameter"
)

# test 17
expect_equal(
  current = sc_object[,, use_cells_to_keep = FALSE],
  target = counts_filtered,
  info = "counts after setting cells to keep and NOT using the parameter"
)

## hvg selection ---------------------------------------------------------------

hvg_to_keep <- 30L

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

# test 18
expect_equivalent(
  current = var_data$mean,
  target = mean_values_r,
  tolerance = 10e-7,
  info = "Correct mean calculations for genes"
)

# test 19
expect_equivalent(
  current = var_data$var,
  target = var_values_r,
  tolerance = 10e-7,
  info = "Correct variance calculations for genes"
)

# due to differences in the loess implementation, this is set higher...
# test 20
expect_equivalent(
  current = var_data$var_std,
  target = var_std,
  tolerance = 10e-3,
  info = "Correct standardised variance calculations for genes"
)

hvg_rs <- get_hvg(sc_object)

# test 21
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

# test 22
expect_equivalent(
  current = var_data$mean,
  target = mean_values_r,
  tolerance = 10e-7,
  info = "Correct mean calculations for genes"
)

# test 23
expect_equivalent(
  current = var_data$var,
  target = var_values_r,
  tolerance = 10e-7,
  info = "Correct variance calculations for genes"
)

# due to differences in the loess implementation, this is set higher...
# test 24
expect_equivalent(
  current = var_data$var_std,
  target = var_std,
  tolerance = 10e-3,
  info = "Correct standardised variance calculations for genes"
)

hvg_rs <- get_hvg(sc_object)

# test 25
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
  no_pcs = 10L,
  randomised_svd = FALSE,
  .verbose = FALSE
)

# test 26
expect_true(
  current = all.equal(
    abs(diag(cor(get_pca_factors(sc_object)[, 1:10], pca_r$x[, 1:10]))),
    rep(1, 10),
    tolerance = 1e-8
  ),
  info = "PCA on single cell data compared to R"
)

sc_object <- calculate_pca_sc(
  object = sc_object,
  no_pcs = 10L,
  randomised_svd = TRUE,
  .verbose = FALSE
)

# test 27
expect_true(
  current = all.equal(
    abs(diag(cor(get_pca_factors(sc_object)[, 1:10], pca_r$x[, 1:10]))),
    rep(1, 10),
    tolerance = 1e-8
  ),
  info = "PCA on single cell data compared to R (randomised SVD)"
)

## knn and snn -----------------------------------------------------------------

if (
  requireNamespace(c("BiocNeighbors", "bluster", "cluster"), quietly = TRUE)
) {
  # annoy algorithm

  sc_object <- find_neighbours_sc(
    sc_object,
    neighbours_params = params_sc_neighbours(knn_algorithm = "annoy"),
    .verbose = FALSE
  )

  expect_true(
    current = class(get_snn_graph(sc_object)) == "igraph",
    info = "igraph correctly returned"
  )

  bioc_knn <- BiocNeighbors::findKNN(
    X = get_pca_factors(sc_object),
    k = 15L
  )$index

  # annoy is a bit more random compared to the implementations in BiocNeighbors
  expect_true(
    current = (sum(sc_object@sc_cache$knn_matrix + 1 == bioc_knn) /
      (dim(bioc_knn)[1] *
        dim(bioc_knn)[2])) >=
      0.98,
    info = "kNN overlap with BiocNeighbors >= 0.98 - annoy algorithm"
  )

  # hnsw
  sc_object <- find_neighbours_sc(
    sc_object,
    neighbours_params = params_sc_neighbours(knn_algorithm = "hnsw"),
    .verbose = FALSE
  )

  expect_true(
    current = (sum(sc_object@sc_cache$knn_matrix + 1 == bioc_knn) /
      (dim(bioc_knn)[1] *
        dim(bioc_knn)[2])) >=
      0.98,
    info = "kNN overlap with BiocNeighbors >= 0.98 - hnsw"
  )

  # snn generation
  bluster_snn <- bluster:::build_snn_graph(t(bioc_knn), "rank", num_threads = 1)

  bixverse_snn <- rs_sc_snn(
    sc_object@sc_cache$knn_matrix,
    snn_method = "rank",
    limited_graph = FALSE,
    pruning = 0,
    verbose = FALSE
  )

  expect_equal(
    current = bixverse_snn$edges + 1,
    target = bluster_snn$edges,
    info = "sNN full generation (rank) between bluster and bixverse - edges"
  )

  expect_true(
    current = cor(bixverse_snn$weights, bluster_snn$weights) >= 0.99,
    info = "sNN full generation (rank) between bluster and bixverse - weights"
  )

  bluster_snn <- bluster:::build_snn_graph(
    t(bioc_knn),
    "jaccard",
    num_threads = 1
  )

  bixverse_snn <- rs_sc_snn(
    sc_object@sc_cache$knn_matrix,
    snn_method = "jaccard",
    limited_graph = FALSE,
    pruning = 0,
    verbose = FALSE
  )

  expect_equal(
    current = bixverse_snn$edges + 1,
    target = bluster_snn$edges,
    info = "sNN full generation (jaccard) between bluster and bixverse - edges"
  )

  expect_true(
    current = cor(bixverse_snn$weights, bluster_snn$weights) >= 0.99,
    info = "sNN full generation (jaccard) between bluster and bixverse - weights"
  )
}

## community detection ---------------------------------------------------------

sc_object <- find_clusters_sc(sc_object)

cell_grps <- unlist(sc_object[["obs_cell_grp"]])
leiden_clusters <- unlist(sc_object[["leiden_clustering"]])

cm <- table(cell_grps, leiden_clusters)

best_match <- apply(cm, 1, which.max)

f1_scores <- sapply(seq_len(nrow(cm)), function(i) {
  tp <- cm[i, best_match[i]]
  fp <- sum(cm[, best_match[i]]) - tp
  fn <- sum(cm[i, ]) - tp

  precision <- tp / (tp + fp)
  recall <- tp / (tp + fn)

  2 * precision * recall / (precision + recall)
})

names(f1_scores) <- rownames(cm)

expect_true(
  current = all(f1_scores > 0.95),
  info = "leiden clustering identifies the cell groups"
)
