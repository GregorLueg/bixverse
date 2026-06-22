# symphony tests ---------------------------------------------------------------

library(magrittr)

test_temp_dir <- file.path(tempdir(), "symphony")
dir.create(test_temp_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot("Test directory does not exist" = dir.exists(test_temp_dir))

## testing parameters ----------------------------------------------------------

# thresholds
min_lib_size <- 300L
min_genes_exp <- 45L
min_cells_exp <- 500L
# hvg / pca / harmony
hvg_to_keep <- 30L
no_pcs <- 15L
harmony_k <- 10L

## helpers ---------------------------------------------------------------------

get_obj_dir <- function(name, current_tmp_dir = test_temp_dir) {
  checkmate::qassert(name, "S1")
  checkmate::assertDirectoryExists(current_tmp_dir)
  final_dir <- file.path(current_tmp_dir, name)
  dir.create(final_dir, showWarnings = FALSE, recursive = TRUE)
  final_dir
}

## synthetic test data ---------------------------------------------------------

# larger data set WITH batch effect
ref_data <- generate_single_cell_test_data(
  syn_data_params = params_sc_synthetic_data(
    n_cells = 1000L,
    n_genes = 100L,
    n_batches = 2L,
    batch_effect_strength = "medium"
  ),
  seed = 1L
)

# smaller data set
query_data <- generate_single_cell_test_data(
  syn_data_params = params_sc_synthetic_data(
    n_cells = 600L,
    n_genes = 100L,
    n_batches = 1L,
    batch_effect_strength = "medium"
  ),
  seed = 2L
)

sc_qc_param_ref <- params_sc_min_quality(
  min_unique_genes = min_genes_exp,
  min_lib_size = min_lib_size,
  min_cells = min_cells_exp,
  target_size = 1000
)

sc_qc_param_query <- params_sc_min_quality(
  min_unique_genes = min_genes_exp,
  min_lib_size = min_lib_size,
  min_cells = 250L,
  target_size = 1000
)

## load reference --------------------------------------------------------------

path_ref <- get_obj_dir(name = "ref")
sc_ref <- SingleCells(dir_data = path_ref)

sc_ref <- load_r_data(
  object = sc_ref,
  counts = ref_data$counts,
  obs = ref_data$obs,
  var = ref_data$var,
  sc_qc_param = sc_qc_param_ref,
  streaming = 0L,
  .verbose = FALSE
)

batch_aware_hvg <- find_hvg_batch_aware_sc(
  object = sc_ref,
  batch_column = "batch_index",
  hvg_no = hvg_to_keep,
  .verbose = FALSE
)

sc_ref <- set_hvg(sc_ref, hvg = batch_aware_hvg$hvg_genes)

ref_hvg <- get_hvg(sc_ref) + 1L # 1-based for the symphony interface

n_hvg_batch_aware <- length(ref_hvg)

## load query ------------------------------------------------------------------

path_query <- get_obj_dir(name = "query")
sc_query <- SingleCells(dir_data = path_query)

sc_query <- load_r_data(
  object = sc_query,
  counts = query_data$counts,
  obs = query_data$obs,
  var = query_data$var,
  sc_qc_param = sc_qc_param_query,
  streaming = 0L,
  .verbose = FALSE
)

# tests ------------------------------------------------------------------------

## build reference (harmony v1) ------------------------------------------------

ref_v1 <- build_symphony_ref(
  object = sc_ref,
  batch_column = "batch_index",
  hvg = ref_hvg,
  harmony_params = params_sc_harmony(k = harmony_k),
  pca_params = params_sc_pca(randomised = FALSE),
  no_pcs = no_pcs,
  seed = 42L,
  .verbose = FALSE
)

n_ref_cells <- length(get_cells_to_keep(sc_ref))

expect_true(
  current = S7::S7_inherits(ref_v1, SymphonyReference),
  info = "build_symphony_ref - returns a SymphonyReference"
)

expect_equal(
  current = S7::prop(ref_v1, "harmony_backend"),
  target = "v1",
  info = "build_symphony_ref - harmony_backend tagged v1"
)

expect_equal(
  current = length(get_symphony_hvg_names(ref_v1)),
  target = n_hvg_batch_aware,
  info = "build_symphony_ref - hvg_gene_names length matches HVG count"
)

expect_equal(
  current = dim(get_symphony_loadings(ref_v1)),
  target = c(n_hvg_batch_aware, no_pcs),
  info = "build_symphony_ref - loadings dimensions"
)

expect_equal(
  current = dim(get_symphony_z_corr(ref_v1)),
  target = c(n_ref_cells, no_pcs),
  info = "build_symphony_ref - z_corr dimensions"
)

expect_equal(
  current = length(S7::prop(ref_v1, "gene_means")),
  target = n_hvg_batch_aware,
  info = "build_symphony_ref - gene_means length"
)

expect_equal(
  current = length(S7::prop(ref_v1, "gene_sds")),
  target = n_hvg_batch_aware,
  info = "build_symphony_ref - gene_sds length"
)

expect_equal(
  current = length(S7::prop(ref_v1, "nr")),
  target = harmony_k,
  info = "build_symphony_ref - nr length matches harmony k"
)

expect_equal(
  current = dim(S7::prop(ref_v1, "centroids")),
  target = c(harmony_k, no_pcs),
  info = "build_symphony_ref - centroids dimensions"
)

## build reference (harmony v2) ------------------------------------------------

ref_v2 <- build_symphony_ref(
  object = sc_ref,
  batch_column = "batch_index",
  hvg = ref_hvg,
  harmony_params = params_sc_harmony_v2(k = harmony_k),
  pca_params = params_sc_pca(randomised = FALSE),
  no_pcs = no_pcs,
  seed = 42L,
  .verbose = FALSE
)

expect_equal(
  current = S7::prop(ref_v2, "harmony_backend"),
  target = "v2",
  info = "build_symphony_ref - harmony_backend tagged v2"
)

expect_equal(
  current = dim(get_symphony_z_corr(ref_v2)),
  target = c(n_ref_cells, no_pcs),
  info = "build_symphony_ref v2 - z_corr dimensions"
)

## slim reference --------------------------------------------------------------

ref_slim <- build_symphony_ref(
  object = sc_ref,
  batch_column = "batch_index",
  hvg = ref_hvg,
  harmony_params = params_sc_harmony(k = harmony_k),
  pca_params = params_sc_pca(randomised = FALSE),
  no_pcs = no_pcs,
  slim = TRUE,
  seed = 42L,
  .verbose = FALSE
)

expect_true(
  current = is.null(S7::prop(ref_slim, "z_orig")),
  info = "slim - z_orig dropped"
)

expect_true(
  current = is.null(S7::prop(ref_slim, "r")),
  info = "slim - r dropped"
)

expect_false(
  current = is.null(S7::prop(ref_slim, "z_corr")),
  info = "slim - z_corr retained for label transfer"
)

expect_true(
  current = S7::prop(ref_slim, "slim"),
  info = "slim - slim flag set"
)

## error: bad harmony_params class ---------------------------------------------

expect_error(
  current = build_symphony_ref(
    object = sc_ref,
    batch_column = "batch_index",
    hvg = ref_hvg,
    harmony_params = list(k = 10), # plain list, no class tag
    pca_params = params_sc_pca(randomised = FALSE),
    no_pcs = no_pcs,
    seed = 42L,
    .verbose = FALSE
  ),
  info = "build_symphony_ref - non-classed harmony_params errors"
)

## map query (with batch) ------------------------------------------------------

sc_query_mapped <- map_symphony_query(
  reference = ref_v1,
  query = sc_query,
  batch_column = "batch_index",
  params = params_symphony_map(),
  .verbose = FALSE
)

n_query_cells <- length(get_cells_to_keep(sc_query_mapped))

z_corr_q <- get_embedding(sc_query_mapped, "symphony")
z_pca_q <- get_embedding(sc_query_mapped, "symphony_pca")
r_q <- get_embedding(sc_query_mapped, "symphony_r")

expect_equal(
  current = dim(z_corr_q),
  target = c(n_query_cells, no_pcs),
  info = "map_symphony_query - z_corr dimensions"
)

expect_equal(
  current = dim(z_pca_q),
  target = c(n_query_cells, no_pcs),
  info = "map_symphony_query - z_pca dimensions"
)

expect_equal(
  current = dim(r_q),
  target = c(n_query_cells, harmony_k),
  info = "map_symphony_query - r dimensions"
)

expect_false(
  current = isTRUE(all.equal(z_corr_q, z_pca_q, check.attributes = FALSE)),
  info = "map_symphony_query - batch correction shifts z_corr away from z_pca"
)

## map query (no batch) --------------------------------------------------------

sc_query_mapped_nb <- map_symphony_query(
  reference = ref_v1,
  query = sc_query,
  batch_column = NULL,
  .verbose = FALSE
)

expect_equivalent(
  current = get_embedding(sc_query_mapped_nb, "symphony"),
  target = get_embedding(sc_query_mapped_nb, "symphony_pca"),
  info = "map_symphony_query - no batch column means z_corr equals z_pca"
)

## label transfer recovers ground truth ----------------------------------------

labels <- transfer_labels_symphony(
  reference = ref_v1,
  reference_object = sc_ref,
  query = sc_query_mapped,
  label_column = "cell_grp",
  knn_params = params_sc_knn(),
  seed = 42L,
  .verbose = FALSE
)

expect_true(
  current = all(
    c("predicted_cell_grp", "confidence_cell_grp") %in% names(labels)
  ),
  info = "transfer_labels_symphony - expected columns present"
)

expect_equal(
  current = nrow(labels),
  target = n_query_cells,
  info = "transfer_labels_symphony - one row per query cell"
)

true_labels <- unlist(sc_query_mapped[["cell_grp"]])
pred_labels <- labels$predicted_cell_grp

accuracy <- mean(true_labels == pred_labels)

expect_true(
  current = accuracy > 0.9,
  info = sprintf(
    "transfer_labels_symphony - accuracy %.1f%% on synthetic data (>90%% expected)",
    100 * accuracy
  )
)

expect_true(
  current = all(
    labels$confidence_cell_grp >= 0 & labels$confidence_cell_grp <= 1
  ),
  info = "transfer_labels_symphony - confidence in [0, 1]"
)

## label transfer with slim reference ------------------------------------------

labels_slim <- transfer_labels_symphony(
  reference = ref_slim,
  reference_object = sc_ref,
  query = sc_query_mapped,
  label_column = "cell_grp",
  knn_params = params_sc_knn(),
  seed = 42L,
  .verbose = FALSE
)

accuracy_slim <- mean(true_labels == labels_slim$predicted_cell_grp)

expect_true(
  current = accuracy_slim > 0.9,
  info = "transfer_labels_symphony - slim reference still recovers labels"
)

# clean up ---------------------------------------------------------------------

on.exit(unlink(test_temp_dir, recursive = TRUE, force = TRUE), add = TRUE)
