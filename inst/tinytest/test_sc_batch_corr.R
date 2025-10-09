# sc batch correction ----------------------------------------------------------

## testing parameters ----------------------------------------------------------

# thresholds
min_lib_size <- 300L
min_genes_exp <- 45L
min_cells_exp <- 1100L
# hvg
hvg_to_keep <- 30L
# pca
no_pcs <- 10L

# slighty more complex data
cell_markers <- list(
  cell_type_1 = list(marker_genes = 0:8L),
  cell_type_2 = list(marker_genes = 9:19L),
  cell_type_3 = list(marker_genes = 20:29L),
  cell_type_4 = list(marker_genes = 30:44L)
)

## synthetic test data ---------------------------------------------------------

# weak batch effect in the data
single_cell_test_data.weak_batch_effect <- generate_single_cell_test_data(
  syn_data_params = params_sc_synthetic_data(
    n_cells = 2000L,
    marker_genes = cell_markers,
    n_batches = 2L,
    batch_effect_strength = "weak"
  )
)

# medium batch effect in the data
single_cell_test_data.medium_batch_effect <- generate_single_cell_test_data(
  syn_data_params = params_sc_synthetic_data(
    n_cells = 2000L,
    marker_genes = cell_markers,
    n_batches = 2L,
    batch_effect_strength = "medium"
  )
)

# strong batch effect in the data
single_cell_test_data.strong_batch_effect <- generate_single_cell_test_data(
  syn_data_params = params_sc_synthetic_data(
    n_cells = 2000L,
    marker_genes = cell_markers,
    n_batches = 2L,
    batch_effect_strength = "strong"
  )
)

# generate QC parameters
sc_qc_param = params_sc_min_quality(
  min_unique_genes = min_genes_exp,
  min_lib_size = min_lib_size,
  min_cells = min_cells_exp
)

# tests ------------------------------------------------------------------------

## pre-processing --------------------------------------------------------------

dir_data <- list("weak" = NULL, "medium" = NULL, "strong" = NULL)
dir_data <- purrr::imap(dir_data, \(elem, name) {
  final_path <- file.path(tempdir(), name)
  dir.create(final_path, showWarnings = FALSE)
  final_path
})

### weak batch effect ----------------------------------------------------------

sc_object.weak_batch_effect <- suppressWarnings(single_cell_exp(
  dir_data = dir_data$weak
))

sc_object.weak_batch_effect <- load_r_data(
  object = sc_object.weak_batch_effect,
  counts = single_cell_test_data.weak_batch_effect$counts,
  obs = single_cell_test_data.weak_batch_effect$obs,
  var = single_cell_test_data.weak_batch_effect$var,
  sc_qc_param = params_sc_min_quality(
    min_unique_genes = min_genes_exp,
    min_lib_size = min_lib_size,
    min_cells = min_cells_exp
  ),
  streaming = FALSE,
  .verbose = FALSE
)

sc_object.weak_batch_effect <- set_cell_to_keep(
  sc_object.weak_batch_effect,
  get_cell_names(sc_object.weak_batch_effect)
)

sc_object.weak_batch_effect <- find_hvg_sc(
  object = sc_object.weak_batch_effect,
  hvg_no = hvg_to_keep,
  .verbose = FALSE
)

sc_object.weak_batch_effect <- calculate_pca_sc(
  object = sc_object.weak_batch_effect,
  no_pcs = no_pcs,
  randomised_svd = FALSE,
  .verbose = FALSE
)

sc_object.weak_batch_effect <- find_neighbours_sc(
  sc_object.weak_batch_effect,
  .verbose = FALSE
)

### medium batch effect --------------------------------------------------------

sc_object.medium_batch_effect <- suppressWarnings(single_cell_exp(
  dir_data = dir_data$medium
))

sc_object.medium_batch_effect <- load_r_data(
  object = sc_object.medium_batch_effect,
  counts = single_cell_test_data.medium_batch_effect$counts,
  obs = single_cell_test_data.medium_batch_effect$obs,
  var = single_cell_test_data.medium_batch_effect$var,
  sc_qc_param = params_sc_min_quality(
    min_unique_genes = min_genes_exp,
    min_lib_size = min_lib_size,
    min_cells = min_cells_exp
  ),
  streaming = FALSE,
  .verbose = FALSE
)

sc_object.medium_batch_effect <- set_cell_to_keep(
  sc_object.medium_batch_effect,
  get_cell_names(sc_object.medium_batch_effect)
)

sc_object.medium_batch_effect <- find_hvg_sc(
  object = sc_object.medium_batch_effect,
  hvg_no = hvg_to_keep,
  .verbose = FALSE
)

sc_object.medium_batch_effect <- calculate_pca_sc(
  object = sc_object.medium_batch_effect,
  no_pcs = no_pcs,
  randomised_svd = FALSE,
  .verbose = FALSE
)

sc_object.medium_batch_effect <- find_neighbours_sc(
  sc_object.medium_batch_effect,
  .verbose = FALSE
)

### strong batch effect --------------------------------------------------------

sc_object.strong_batch_effect <- suppressWarnings(single_cell_exp(
  dir_data = dir_data$strong
))

sc_object.strong_batch_effect <- load_r_data(
  object = sc_object.strong_batch_effect,
  counts = single_cell_test_data.strong_batch_effect$counts,
  obs = single_cell_test_data.strong_batch_effect$obs,
  var = single_cell_test_data.strong_batch_effect$var,
  sc_qc_param = params_sc_min_quality(
    min_unique_genes = min_genes_exp,
    min_lib_size = min_lib_size,
    min_cells = min_cells_exp
  ),
  streaming = FALSE,
  .verbose = FALSE
)

sc_object.strong_batch_effect <- set_cell_to_keep(
  sc_object.strong_batch_effect,
  get_cell_names(sc_object.strong_batch_effect)
)

sc_object.strong_batch_effect <- find_hvg_sc(
  object = sc_object.strong_batch_effect,
  hvg_no = hvg_to_keep,
  .verbose = FALSE
)

sc_object.strong_batch_effect <- calculate_pca_sc(
  object = sc_object.strong_batch_effect,
  no_pcs = no_pcs,
  randomised_svd = FALSE,
  .verbose = FALSE
)

sc_object.strong_batch_effect <- find_neighbours_sc(
  sc_object.strong_batch_effect,
  .verbose = FALSE
)

# tests ------------------------------------------------------------------------

## kbet scores -----------------------------------------------------------------

kbet_scores.weak_batch_effect <- calculate_kbet_sc(
  object = sc_object.weak_batch_effect,
  batch_column = "batch_index"
)

kbet_scores.medium_batch_effect <- calculate_kbet_sc(
  object = sc_object.medium_batch_effect,
  batch_column = "batch_index"
)

kbet_scores.strong_batch_effect <- calculate_kbet_sc(
  object = sc_object.strong_batch_effect,
  batch_column = "batch_index"
)

expect_true(
  current = is.logical(kbet_scores.weak_batch_effect$significant_tests),
  info = paste("kbet scores - returns booleans where expected")
)

expect_true(
  current = checkmate::qtest(
    kbet_scores.weak_batch_effect$chisquare_pvals,
    "N[0, 1]"
  ),
  info = paste("kbet scores - p-values are expected type and range")
)

expect_true(
  current = kbet_scores.weak_batch_effect$kbet_score <
    kbet_scores.medium_batch_effect$kbet_score,
  info = paste("kbet scores - weak kbet < medium kbet")
)

expect_true(
  current = kbet_scores.weak_batch_effect$kbet_score <
    kbet_scores.strong_batch_effect$kbet_score,
  info = paste("kbet scores - weak kbet < strong kbet")
)

expect_true(
  current = kbet_scores.medium_batch_effect$kbet_score <
    kbet_scores.strong_batch_effect$kbet_score,
  info = paste("kbet scores - weak kbet < strong kbet")
)
