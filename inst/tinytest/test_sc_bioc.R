# comparisons with blusters, BiocNeighbors -------------------------------------

if (
  !all(sapply(
    c("BiocNeighbors", "bluster", "cluster"),
    requireNamespace,
    quietly = TRUE
  ))
) {
  exit_file("BiocNeighbors, bluster and cluster not available")
}

## test parameters -------------------------------------------------------------

# thresholds
min_lib_size <- 300L
min_genes_exp <- 45L
min_cells_exp <- 500L
hvg_to_keep <- 30L
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


## quick run -------------------------------------------------------------------

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

sc_qc_param = params_sc_min_quality(
  min_unique_genes = min_genes_exp,
  min_lib_size = min_lib_size,
  min_cells = min_cells_exp,
  target_size = 1000
)

sc_object <- suppressWarnings(single_cell_exp(dir_data = tempdir()))

sc_object <- load_h5ad(
  object = sc_object,
  h5_path = path.expand(f_path_csr),
  sc_qc_param = sc_qc_param,
  .verbose = FALSE
)

gs_of_interest <- list(
  gs_1 = c("gene_001", "gene_002", "gene_003", "gene_004"),
  gs_2 = c("gene_096", "gene_097", "gene_100")
)

sc_object <- gene_set_proportions_sc(
  sc_object,
  gs_of_interest,
  .verbose = FALSE
)

threshold <- 0.05

cells_to_keep <- sc_object[[]][gs_2 < threshold, cell_id]

expect_true(
  current = length(cells_to_keep) > 600,
  info = "sensible cell filtering based on the threshold"
)

sc_object <- set_cell_to_keep(sc_object, cells_to_keep)

sc_object <- find_hvg_sc(
  object = sc_object,
  hvg_no = hvg_to_keep,
  .verbose = FALSE
)

sc_object <- calculate_pca_sc(
  object = sc_object,
  no_pcs = no_pcs,
  randomised_svd = FALSE,
  .verbose = FALSE
)

# annoy algorithm

sc_object <- find_neighbours_sc(
  sc_object,
  neighbours_params = params_sc_neighbours(knn_algorithm = "annoy"),
  .verbose = FALSE
)

bioc_knn <- BiocNeighbors::findKNN(
  X = get_pca_factors(sc_object),
  k = 15L
)$index


# annoy is a bit more random compared to the implementations in BiocNeighbors
expect_true(
  current = (sum(get_knn_mat(sc_object) + 1 == bioc_knn) /
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
  current = (sum(get_knn_mat(sc_object) + 1 == bioc_knn) /
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
