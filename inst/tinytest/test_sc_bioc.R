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

test_temp_dir <- file.path(
  tempdir(),
  paste0("test_", format(Sys.time(), "%Y%m%d_%H%M%S_"), sample(1000:9999, 1))
)
dir.create(test_temp_dir, recursive = TRUE)

## test parameters -------------------------------------------------------------

# thresholds
min_lib_size <- 300L
min_genes_exp <- 45L
min_cells_exp <- 500L
hvg_to_keep <- 50L
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

sc_qc_param = params_sc_min_quality(
  min_unique_genes = min_genes_exp,
  min_lib_size = min_lib_size,
  min_cells = min_cells_exp,
  target_size = 1000
)

sc_object <- single_cell_exp(dir_data = test_temp_dir)

sc_object <- load_r_data(
  object = sc_object,
  counts = single_cell_test_data$counts,
  obs = single_cell_test_data$obs,
  var = single_cell_test_data$var,
  sc_qc_param = sc_qc_param,
  .verbose = FALSE
)

sc_object <- set_cells_to_keep(sc_object, get_cell_names(sc_object))

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

#### bioc ----------------------------------------------------------------------

bioc_knn_cosine <- BiocNeighbors::findKNN(
  X = get_pca_factors(sc_object),
  k = 15L,
  BNPARAM = BiocNeighbors::AnnoyParam(distance = "Cosine")
)$index

bioc_knn_euclidean <- BiocNeighbors::findKNN(
  X = get_pca_factors(sc_object),
  k = 15L
)$index

#### annoy cosine --------------------------------------------------------------

# cosine
sc_object <- find_neighbours_sc(
  sc_object,
  neighbours_params = params_sc_neighbours(
    knn = list(k = 15L, knn_method = "annoy")
  ),
  .verbose = FALSE
)

sc_knn_cosine <- get_knn_mat(sc_object)

# annoy is a bit more random compared to the implementations in BiocNeighbors
expect_true(
  current = (sum(sc_knn_cosine + 1 == bioc_knn_cosine) /
    (prod(dim(bioc_knn_cosine)))) >=
    0.95,
  info = paste(
    "kNN overlap with BiocNeighbors >= 0.95",
    "- annoy algorithm (cosine)"
  )
)

#### annoy euclidean -----------------------------------------------------------

sc_object <- find_neighbours_sc(
  sc_object,
  neighbours_params = params_sc_neighbours(
    knn = list(knn_method = "annoy", ann_dist = "euclidean", k = 15L)
  ),
  .verbose = FALSE
)

sc_knn_euclidean <- get_knn_mat(sc_object)

expect_true(
  current = (sum(sc_knn_euclidean + 1 == bioc_knn_euclidean) /
    (prod(dim(bioc_knn_euclidean)))) >=
    0.95,
  info = paste(
    "kNN overlap with BiocNeighbors >= 0.95",
    "- annoy algorithm (euclidean)"
  )
)

#### hnsw cosine ---------------------------------------------------------------

sc_object <- find_neighbours_sc(
  sc_object,
  neighbours_params = params_sc_neighbours(
    knn = list(knn_method = "hnsw", k = 15L)
  ),
  .verbose = FALSE
)

sc_knn_cosine <- get_knn_mat(sc_object)

expect_true(
  current = (sum(sc_knn_cosine + 1 == bioc_knn_cosine) /
    (prod(dim(bioc_knn_cosine)))) >=
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
    knn = list(
      knn_algorithm = "hnsw",
      ann_dist = "euclidean",
      k = 15L
    )
  ),
  .verbose = FALSE
)

sc_knn_euclidean <- get_knn_mat(sc_object)

expect_true(
  current = (sum(sc_knn_euclidean + 1 == bioc_knn_euclidean) /
    (prod(dim(bioc_knn_euclidean)))) >=
    0.95,
  info = paste(
    "kNN overlap with BiocNeighbors >= 0.9",
    "- HNSW algorithm (euclidean)"
  )
)

#### nndescent cosine ----------------------------------------------------------

# NNDescent overlaps way worse than the others due to the underlying
# synthetic data... On the full data it should still reach 0.75 and on TopX
# neighbours better...

sc_object <- find_neighbours_sc(
  sc_object,
  neighbours_params = params_sc_neighbours(
    knn = list(
      knn_method = "nndescent",
      k = 15L
    )
  ),
  .verbose = FALSE
)

sc_knn_cosine <- get_knn_mat(sc_object)

nn_descent_full_overlap <- (sum(sc_knn_cosine + 1 == bioc_knn_cosine) /
  (dim(bioc_knn_cosine)[1] *
    dim(bioc_knn_cosine)[2]))

nn_descent_top5_overlap <- (sum(
  (sc_knn_cosine + 1)[, 1:5] == bioc_knn_cosine[, 1:5]
) /
  (dim(bioc_knn_cosine)[1] * 5))

expect_true(
  current = nn_descent_full_overlap >= 0.75,
  info = paste(
    "kNN overlap with BiocNeighbors >= 0.75",
    "- NNDescent algorithm (cosine)"
  )
)

expect_true(
  current = nn_descent_top5_overlap >= 0.95,
  info = paste(
    "kNN overlap with BiocNeighbors >= 0.95",
    "- NNDescent algorithm (cosine) Top 5 neighbours"
  )
)

#### nndescent euclidean -------------------------------------------------------

sc_object <- find_neighbours_sc(
  sc_object,
  neighbours_params = params_sc_neighbours(
    knn = list(
      knn_algorithm = "nndescent",
      ann_dist = "euclidean",
      k = 15L
    )
  ),
  .verbose = FALSE
)

sc_knn_euclidean <- get_knn_mat(sc_object)

nn_descent_full_overlap <- (sum(sc_knn_euclidean + 1 == bioc_knn_euclidean) /
  (dim(bioc_knn_euclidean)[1] *
    dim(bioc_knn_euclidean)[2]))

nn_descent_top5_overlap <- (sum(
  (sc_knn_euclidean + 1)[, 1:5] == bioc_knn_euclidean[, 1:5]
) /
  (dim(bioc_knn_euclidean)[1] * 5))

expect_true(
  current = nn_descent_full_overlap >= 0.75,
  info = paste(
    "kNN overlap with BiocNeighbors >= 0.75",
    "- NNDescent algorithm (euclidean)"
  )
)

expect_true(
  current = nn_descent_top5_overlap >= 0.95,
  info = paste(
    "kNN overlap with BiocNeighbors >= 0.75",
    "- NNDescent algorithm (euclidean) Top 5 neighbours"
  )
)

### snn graph ------------------------------------------------------------------

#### bluster -------------------------------------------------------------------

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

#### hnsw euclidean ------------------------------------------------------------

# due to tiny differences, i will only test the case where the similarity
# is perfect, i.e., HNSW == 1.

sc_object <- find_neighbours_sc(
  sc_object,
  neighbours_params = params_sc_neighbours(
    knn = list(
      knn_method = "hnsw",
      ann_dist = "euclidean"
    )
  ),
  .verbose = FALSE
)

sc_knn_euclidean <- get_knn_mat(sc_object)

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

## generate knn with distances -------------------------------------------------

cells_to_keep <- sc_object[[]][cell_grp != "cell_type_1", cell_id]

### class and getters ----------------------------------------------------------

annoy_knn_data <- generate_knn_sc(
  sc_object,
  neighbours_params = params_sc_neighbours(
    knn = list(
      knn_method = "annoy",
      ann_dist = "euclidean"
    )
  ),
  .verbose = FALSE
)

expect_true(
  current = checkmate::testClass(annoy_knn_data, "sc_knn"),
  info = "sc_knn - expected class returned"
)

expect_true(
  current = checkmate::testMatrix(
    get_knn_mat(annoy_knn_data),
    mode = "integer"
  ),
  info = "sc_knn - indices of expected type"
)

expect_true(
  current = checkmate::testMatrix(
    get_knn_dist(annoy_knn_data),
    mode = "numeric"
  ),
  info = "sc_knn - distances of expected type"
)

annoy_knn_data_red <- generate_knn_sc(
  sc_object,
  cells_to_use = cells_to_keep,
  neighbours_params = params_sc_neighbours(
    knn = list(
      knn_method = "annoy",
      ann_dist = "euclidean"
    )
  ),
  .verbose = FALSE
)

expect_equivalent(
  current = annoy_knn_data_red$used_cells,
  target = cells_to_keep,
  info = "sc_knn - reduced version behaves"
)

### annoy ----------------------------------------------------------------------

annoy_knn_data_cosine <- generate_knn_sc(
  sc_object,
  neighbours_params = params_sc_neighbours(
    knn = list(
      knn_method = "annoy",
      ann_dist = "cosine"
    )
  ),
  .verbose = FALSE
)

sc_knn_annoy_euclidean_val <- sum(
  (get_knn_mat(annoy_knn_data) + 1) == bioc_knn_euclidean
) /
  prod(dim(bioc_knn_euclidean))

sc_knn_annoy_cosine_val <- sum(
  (get_knn_mat(annoy_knn_data_cosine) + 1) == bioc_knn_cosine
) /
  prod(dim(bioc_knn_euclidean))

expect_true(
  current = sc_knn_annoy_euclidean_val >= 0.95,
  info = "sc_knn class - expected overlap annoy euclidean"
)

expect_true(
  current = sc_knn_annoy_cosine_val >= 0.95,
  info = "sc_knn class - expected overlap annoy cosine"
)

### hnsw -----------------------------------------------------------------------

hnsw_knn_data <- generate_knn_sc(
  sc_object,
  neighbours_params = params_sc_neighbours(
    knn = list(
      knn_method = "hnsw",
      ann_dist = "euclidean"
    )
  ),
  .verbose = FALSE
)

sc_knn_hnsw_euclidean_val <- sum(
  (get_knn_mat(hnsw_knn_data) + 1) == bioc_knn_euclidean
) /
  prod(dim(bioc_knn_euclidean))

expect_true(
  current = sc_knn_hnsw_euclidean_val >= 0.95,
  info = "sc_knn class - expected overlap hnsw euclidean"
)

### nndescent ------------------------------------------------------------------

nndescent_knn_data <- generate_knn_sc(
  sc_object,
  neighbours_params = params_sc_neighbours(
    knn = list(
      knn_method = "nndescent",
      ann_dist = "euclidean"
    )
  ),
  .verbose = FALSE
)

sc_knn_nndescent_euclidean_val <- sum(
  (get_knn_mat(nndescent_knn_data) + 1) == bioc_knn_euclidean
) /
  prod(dim(bioc_knn_euclidean))

sc_nndescent_top5_overlap <- (sum(
  (get_knn_mat(nndescent_knn_data) + 1)[, 1:5] == bioc_knn_euclidean[, 1:5]
) /
  (dim(bioc_knn_euclidean)[1] * 5))

# worse overlap here... as expected
expect_true(
  current = sc_knn_nndescent_euclidean_val >= 0.75,
  info = "sc_knn class - expected overlap nndescent euclidean"
)

# worse overlap here... as expected
expect_true(
  current = sc_nndescent_top5_overlap >= 0.95,
  info = "sc_knn class - expected overlap nndescent euclidean - top5 neighbours"
)
