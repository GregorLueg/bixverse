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
no_pcs <- 15L

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

## direct rust <> bioc comparison ----------------------------------------------

pca_embd <- get_pca_factors(sc_object)

# Helper function to compare against Bioc results
# Needs to take care of 1-indexing vs 0-indexing
calc_recall_bioc <- function(knn_mat, rs_knn_mat) {
  sum(knn_mat == (rs_knn_mat + 1)) / prod(dim(knn_mat))
}

# Helper function to calculate the Jaccard similarity
call_jaccard_bioc <- function(knn_mat, rs_knn_mat) {
  rs_jaccard_row_integers(knn_mat, rs_knn_mat + 1L)
}

# Helper function to calculate against exhaustive search
calc_recall_exh <- function(exhaustive_search_knn, knn_mat) {
  sum(exhaustive_search_knn == knn_mat) / prod(dim(knn_mat))
}

# Helper function to calculate the Jaccard similarity against exhaustive search
call_jaccard_exh <- function(exhaustive_search_knn, rs_knn_mat) {
  rs_jaccard_row_integers(exhaustive_search_knn, rs_knn_mat)
}

### bio part -------------------------------------------------------------------

bioc_knn_euclidean <- BiocNeighbors::findKNN(
  X = get_pca_factors(sc_object),
  k = 15L
)$index

bioc_knn_cosine <- BiocNeighbors::findKNN(
  X = get_pca_factors(sc_object),
  k = 15L,
  BNPARAM = BiocNeighbors::AnnoyParam(distance = "Cosine")
)$index


### rust part ------------------------------------------------------------------

#### euclidean distances -------------------------------------------------------

rs_exhaustive_euc <- rs_sc_knn(
  embd = pca_embd,
  knn_params = list(knn_method = "exhaustive", ann_dist = "euclidean"),
  verbose = FALSE,
  seed = 42L
)

rs_annoy_euc <- rs_sc_knn(
  embd = pca_embd,
  knn_params = list(knn_method = "annoy", ann_dist = "euclidean"),
  verbose = FALSE,
  seed = 42L
)

rs_hnsw_euc <- rs_sc_knn(
  embd = pca_embd,
  knn_params = list(knn_method = "hnsw", ann_dist = "euclidean"),
  verbose = FALSE,
  seed = 42L
)

rs_nndescent_euc <- rs_sc_knn(
  embd = pca_embd,
  knn_params = list(knn_method = "nndescent", ann_dist = "euclidean"),
  verbose = FALSE,
  seed = 42L
)

##### comparison against bioc --------------------------------------------------

bioc_exhaustive_euc <- calc_recall_bioc(bioc_knn_euclidean, rs_exhaustive_euc)
bioc_annoy_euc <- calc_recall_bioc(bioc_knn_euclidean, rs_annoy_euc)
bioc_hnsw_euc <- calc_recall_bioc(bioc_knn_euclidean, rs_hnsw_euc)
bioc_nndescent_euc <- calc_recall_bioc(bioc_knn_euclidean, rs_nndescent_euc)

expect_true(
  current = bioc_exhaustive_euc >= 0.98,
  info = "bioc <> exhaustive search (euclidean dist) - high overlap"
)
expect_true(
  current = bioc_annoy_euc >= 0.98,
  info = "bioc <> Annoy (euclidean dist) - high overlap"
)
expect_true(
  current = bioc_hnsw_euc >= 0.98,
  info = "bioc <> HNSW (euclidean dist) - high overlap"
)
expect_true(
  current = bioc_nndescent_euc >= 0.98,
  info = "bioc <> NNDescent (euclidean dist) - high overlap"
)

##### against ground truth -----------------------------------------------------

gt_annoy_euc <- calc_recall_exh(rs_exhaustive_euc, rs_annoy_euc)
gt_hnsw_euc <- calc_recall_exh(rs_exhaustive_euc, rs_hnsw_euc)
gt_nndescent_euc <- calc_recall_exh(rs_exhaustive_euc, rs_nndescent_euc)

expect_true(
  current = gt_annoy_euc >= 0.98,
  info = "annoy <> exhaustive search (euclidean dist) - high overlap"
)
expect_true(
  current = gt_hnsw_euc >= 0.98,
  info = "hnsw <> exhaustive search (euclidean dist) - high overlap"
)
expect_true(
  current = gt_nndescent_euc >= 0.98,
  info = "nndescent <> exhaustive search (euclidean dist) - high overlap"
)

#### cosine distances ----------------------------------------------------------

rs_exhaustive_cos <- rs_sc_knn(
  embd = pca_embd,
  knn_params = list(knn_method = "exhaustive", ann_dist = "cosine"),
  verbose = FALSE,
  seed = 42L
)

rs_annoy_cos <- rs_sc_knn(
  embd = pca_embd,
  knn_params = list(knn_method = "annoy", ann_dist = "cosine"),
  verbose = FALSE,
  seed = 42L
)

rs_hnsw_cos <- rs_sc_knn(
  embd = pca_embd,
  knn_params = list(knn_method = "hnsw", ann_dist = "cosine"),
  verbose = FALSE,
  seed = 42L
)

rs_nndescent_cos <- rs_sc_knn(
  embd = pca_embd,
  knn_params = list(knn_method = "nndescent", ann_dist = "cosine"),
  verbose = FALSE,
  seed = 42L
)

##### comparison against bioc --------------------------------------------------

# in all cases, there is a slight difference... the exhaustive difference
# however is correct...

bioc_exhaustive_cos <- calc_recall_bioc(bioc_knn_cosine, rs_exhaustive_cos)
bioc_annoy_cos <- calc_recall_bioc(bioc_knn_cosine, rs_annoy_cos)
bioc_hnsw_cos <- calc_recall_bioc(bioc_knn_cosine, rs_hnsw_cos)
bioc_nndescent_cos <- calc_recall_bioc(bioc_knn_cosine, rs_hnsw_cos)

expect_true(
  current = bioc_exhaustive_cos >= 0.98,
  info = "bioc <> exhaustive search (cosine dist) - high overlap"
)
expect_true(
  current = bioc_annoy_cos >= 0.98,
  info = "bioc <> Annoy (cosine dist) - high overlap"
)
expect_true(
  current = bioc_hnsw_cos >= 0.98,
  info = "bioc <> HNSW (cosine dist) - high overlap"
)
expect_true(
  current = bioc_nndescent_cos >= 0.98,
  info = "bioc <> NNDescent (cosine dist) - high overlap"
)

##### against ground truth -----------------------------------------------------

gt_annoy_cos <- calc_recall_exh(rs_exhaustive_cos, rs_annoy_cos)
gt_hnsw_cos <- calc_recall_exh(rs_exhaustive_cos, rs_hnsw_cos)
gt_nndescent_cos <- calc_recall_exh(rs_exhaustive_cos, rs_nndescent_cos)

expect_true(
  current = gt_annoy_cos >= 0.98,
  info = "annoy <> exhaustive search (cosine dist) - high overlap"
)
expect_true(
  current = gt_hnsw_cos >= 0.98,
  info = "hnsw <> exhaustive search (cosine dist) - high overlap"
)
expect_true(
  current = gt_nndescent_cos >= 0.98,
  info = "nndescent <> exhaustive search (cosine dist) - high overlap"
)

## object ----------------------------------------------------------------------

### parameters -----------------------------------------------------------------

# overwrite defaults
non_default_params <- params_sc_neighbours(
  knn = list(k = 10L, knn_method = "annoy", ann_dist = "euclidean")
)

expect_true(
  current = non_default_params$knn_method == "annoy",
  info = "neighbour default params are overwritten - knn method"
)

expect_true(
  current = non_default_params$k == 10L,
  info = "neighbour default params are overwritten - k"
)

expect_true(
  current = non_default_params$ann_dist == "euclidean",
  info = "neighbour default params are overwritten - distance measure"
)

### correct data is returned ---------------------------------------------------

#### check consistency with rust -----------------------------------------------

# default with HNSW and cosine
sc_object <- find_neighbours_sc(
  sc_object,
  neighbours_params = params_sc_neighbours(
    knn = list(k = 15L, knn_method = "hnsw")
  ),
  .verbose = FALSE
)

sc_knn_cosine <- get_knn_mat(sc_object)

expect_equivalent(
  current = sc_knn_cosine,
  target = rs_hnsw_cos,
  info = "method wrapper behaving for HNSW"
)

# default with HNSW and cosine
sc_object <- find_neighbours_sc(
  sc_object,
  neighbours_params = params_sc_neighbours(
    knn = list(k = 15L, knn_method = "annoy", ann_dist = "euclidean")
  ),
  .verbose = FALSE
)

sc_knn_euclidean <- get_knn_mat(sc_object)

expect_equivalent(
  current = sc_knn_euclidean,
  target = rs_annoy_euc,
  info = "method wrapper behaving for Annoy with different distance metric"
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


expect_true(
  current = calc_recall_exh(rs_annoy_euc, get_knn_mat(annoy_knn_data)) == 1,
  info = "sc_knn class - expected overlap annoy euclidean"
)

expect_true(
  current = calc_recall_exh(rs_annoy_cos, get_knn_mat(annoy_knn_data_cosine)) ==
    1,
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

expect_true(
  current = calc_recall_exh(rs_hnsw_euc, get_knn_mat(hnsw_knn_data)) == 1,
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

expect_true(
  current = calc_recall_exh(
    rs_nndescent_euc,
    get_knn_mat(nndescent_knn_data)
  ) ==
    1,
  info = "sc_knn class - expected overlap nndescent euclidean"
)
