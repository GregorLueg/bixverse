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
no_pcs <- 5L

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

### knn generation -------------------------------------------------------------

#### annoy cosine --------------------------------------------------------------

# Bioc results
bioc_knn_cosine <- BiocNeighbors::findKNN(
  X = get_pca_factors(sc_object),
  k = 15L,
  BNPARAM = BiocNeighbors::AnnoyParam(distance = "Cosine")
)$index

# cosine
sc_object <- find_neighbours_sc(
  sc_object,
  neighbours_params = params_sc_neighbours(
    knn_algorithm = "annoy"
  ),
  .verbose = FALSE
)

sc_knn_cosine <- get_knn_mat(sc_object)

# annoy is a bit more random compared to the implementations in BiocNeighbors
expect_true(
  current = (sum(sc_knn_cosine + 1 == bioc_knn_cosine) /
    (dim(bioc_knn_cosine)[1] *
      dim(bioc_knn_cosine)[2])) >=
    0.95,
  info = paste(
    "kNN overlap with BiocNeighbors >= 0.95",
    "- annoy algorithm (cosine)"
  )
)

#### annoy euclidean -----------------------------------------------------------

bioc_knn_euclidean <- BiocNeighbors::findKNN(
  X = get_pca_factors(sc_object),
  k = 15L
)$index

sc_object <- find_neighbours_sc(
  sc_object,
  neighbours_params = params_sc_neighbours(
    knn_algorithm = "annoy",
    ann_dist = "euclidean"
  ),
  .verbose = FALSE
)

sc_knn_euclidean <- get_knn_mat(sc_object)

expect_true(
  current = (sum(sc_knn_euclidean + 1 == bioc_knn_euclidean) /
    (dim(bioc_knn_euclidean)[1] *
      dim(bioc_knn_euclidean)[2])) >=
    0.95,
  info = paste(
    "kNN overlap with BiocNeighbors >= 0.95",
    "- annoy algorithm (euclidean)"
  )
)

#### hnsw cosine ---------------------------------------------------------------

sc_object <- find_neighbours_sc(
  sc_object,
  neighbours_params = params_sc_neighbours(knn_algorithm = "hnsw"),
  .verbose = FALSE
)

sc_knn_cosine <- get_knn_mat(sc_object)

expect_true(
  current = (sum(sc_knn_cosine + 1 == bioc_knn_cosine) /
    (dim(bioc_knn_cosine)[1] *
      dim(bioc_knn_cosine)[2])) >=
    0.95,
  info = paste(
    "kNN overlap with BiocNeighbors >= 0.95",
    "- HNSW algorithm (cosine)"
  )
)

#### hnsw euclidean ------------------------------------------------------------

sc_object <- find_neighbours_sc(
  sc_object,
  neighbours_params = params_sc_neighbours(
    knn_algorithm = "hnsw",
    ann_dist = "euclidean"
  ),
  .verbose = FALSE
)

sc_knn_euclidean <- get_knn_mat(sc_object)

expect_true(
  current = (sum(sc_knn_euclidean + 1 == bioc_knn_euclidean) /
    (dim(bioc_knn_euclidean)[1] *
      dim(bioc_knn_euclidean)[2])) >=
    0.95,
  info = paste(
    "kNN overlap with BiocNeighbors >= 0.95",
    "- HNSW algorithm (euclidean)"
  )
)

### snn graph ------------------------------------------------------------------

#### euclidean -----------------------------------------------------------------

# snn generation
bluster_snn_euclidean_rank <- bluster:::build_snn_graph(
  t(bioc_knn_euclidean),
  "rank",
  num_threads = 1
)

bluster_snn_euclidean_jaccard <- bluster:::build_snn_graph(
  t(bioc_knn_euclidean),
  "jaccard",
  num_threads = 1
)

bixverse_snn_euclidean_rank <- rs_sc_snn(
  sc_knn_euclidean,
  snn_method = "rank",
  limited_graph = FALSE,
  pruning = 0,
  verbose = FALSE
)

bixverse_snn_euclidean_jaccard <- rs_sc_snn(
  sc_knn_euclidean,
  snn_method = "jaccard",
  limited_graph = FALSE,
  pruning = 0,
  verbose = FALSE
)

# rank version
expect_equal(
  current = bixverse_snn_euclidean_rank$edges + 1,
  target = bluster_snn_euclidean_rank$edges,
  info = paste(
    "sNN full generation (euclidean; rank)",
    "between bluster and bixverse - edges"
  )
)

expect_true(
  current = cor(
    bixverse_snn_euclidean_rank$weights,
    bluster_snn_euclidean_rank$weights
  ) >=
    0.95,
  info = paste(
    "sNN full generation (euclidean; rank)",
    "between bluster and bixverse - weights"
  )
)

# jaccard
expect_equal(
  current = bixverse_snn_euclidean_jaccard$edges + 1,
  target = bluster_snn_euclidean_jaccard$edges,
  info = paste(
    "sNN full generation (euclidean; jaccard)",
    "between bluster and bixverse - edges"
  )
)

expect_true(
  current = cor(
    bixverse_snn_euclidean_jaccard$weights,
    bluster_snn_euclidean_jaccard$weights
  ) >=
    0.95,
  info = paste(
    "sNN full generation (euclidean; jaccard)",
    "between bluster and bixverse - weights"
  )
)

# #### cosine --------------------------------------------------------------------

# snn generation
bluster_snn_cosine_rank <- bluster:::build_snn_graph(
  t(bioc_knn_cosine),
  "rank",
  num_threads = 1
)

bluster_snn_cosine_jaccard <- bluster:::build_snn_graph(
  t(bioc_knn_cosine),
  "jaccard",
  num_threads = 1
)

bixverse_snn_cosine_rank <- rs_sc_snn(
  sc_knn_cosine,
  snn_method = "rank",
  limited_graph = FALSE,
  pruning = 0,
  verbose = FALSE
)

bixverse_snn_cosine_jaccard <- rs_sc_snn(
  sc_knn_cosine,
  snn_method = "jaccard",
  limited_graph = FALSE,
  pruning = 0,
  verbose = FALSE
)

# rank version

expect_equal(
  current = bixverse_snn_cosine_rank$edges + 1,
  target = bluster_snn_cosine_rank$edges,
  info = paste(
    "sNN full generation (euclidean; rank)",
    "between bluster and bixverse - edges"
  )
)

expect_true(
  current = cor(
    bixverse_snn_cosine_rank$weights,
    bluster_snn_cosine_rank$weights
  ) >=
    0.99,
  info = paste(
    "sNN full generation (euclidean; rank)",
    "between bluster and bixverse - weights"
  )
)
