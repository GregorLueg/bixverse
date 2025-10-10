# sc batch correction ----------------------------------------------------------

## testing parameters ----------------------------------------------------------

# thresholds
min_lib_size <- 300L
min_genes_exp <- 45L
min_cells_exp <- 750L
# hvg
hvg_to_keep <- 30L
# pca
no_pcs <- 15L

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
    n_cells = 1500L,
    marker_genes = cell_markers,
    n_batches = 3L,
    batch_effect_strength = "weak"
  )
)

# medium batch effect in the data
single_cell_test_data.medium_batch_effect <- generate_single_cell_test_data(
  syn_data_params = params_sc_synthetic_data(
    n_cells = 1500L,
    marker_genes = cell_markers,
    n_batches = 3L,
    batch_effect_strength = "medium"
  )
)

# strong batch effect in the data
single_cell_test_data.strong_batch_effect <- generate_single_cell_test_data(
  syn_data_params = params_sc_synthetic_data(
    n_cells = 1500L,
    marker_genes = cell_markers,
    n_batches = 3L,
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
  neighbours_params = params_sc_neighbours(k = 10L),
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
  neighbours_params = params_sc_neighbours(k = 10L),
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
  neighbours_params = params_sc_neighbours(k = 10L),
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

## bbknn -----------------------------------------------------------------------

### helper function for assessment of batch correction

assess_bbknn_impact <- function(object) {
  # checks
  checkmate::assertTRUE(S7::S7_inherits(object, single_cell_exp))

  # get key data from the object
  batch_labels <- unlist(object[["batch_index"]])
  emdb_mat <- get_pca_factors(object)
  knn_original <- get_knn_mat(object)
  cell_types <- unlist(object[["cell_grp"]])
  batch_effects <- unlist(object[["batch_index"]])

  checkmate::assertTRUE(nrow(emdb_mat) == length(batch_effects))

  # get the original kbet score and see how many cells are neighbour
  # of same cell type as itself
  original_kbet <- kbet_scores.weak_batch_effect <- calculate_kbet_sc(
    object = object,
    batch_column = "batch_index"
  )

  correct_neighbours_uncor <- vector(
    mode = "numeric",
    length = length(cell_types)
  )

  for (cell_idx in seq_len(length(cell_types))) {
    correct_neighbours_uncor[cell_idx] <- sum(
      cell_types[knn_original[cell_idx, ] + 1] == cell_types[cell_idx]
    )
  }

  bbknn_corrected_data <- rs_bbknn(
    embd = emdb_mat,
    batch_labels = as.integer(batch_labels),
    bbknn_params = list(
      knn_method = "annoy",
      neighbours_within_batch = 5L,
      set_op_mix_ratio = 0.5
    ),
    seed = 42L,
    verbose = FALSE
  )

  knn_corr <- rs_bbknn_filtering(
    indptr = bbknn_corrected_data$distances$indptr,
    indices = bbknn_corrected_data$distances$indices,
    10L
  )

  storage.mode(knn_corr) <- "integer"

  k_bet_corrected <- sum(rs_kbet(knn_corr, as.integer(batch_effects)) < 0.05)

  correct_neighbours_cor <- vector(
    mode = "numeric",
    length = length(cell_types)
  )

  for (cell_idx in seq_len(length(cell_types))) {
    correct_neighbours_cor[cell_idx] <- sum(
      cell_types[knn_corr[cell_idx, ] + 1] == cell_types[cell_idx]
    )
  }

  res <- list(
    kbet_original = original_kbet$kbet_score,
    kbet_correct = k_bet_corrected,
    neighbours_original = correct_neighbours_uncor,
    neighbours_corrected = correct_neighbours_cor
  )

  return(res)
}

### general logic of the rust code ---------------------------------------------

#### weak batch effects --------------------------------------------------------

weak_batch_effect_res <- assess_bbknn_impact(
  object = sc_object.weak_batch_effect
)

expect_true(
  current = weak_batch_effect_res$kbet_original >
    weak_batch_effect_res$kbet_correct,
  info = paste("bbknn - weak batch effects get regressed out")
)

expect_true(
  current = mean(weak_batch_effect_res$neighbours_corrected) >
    (mean(weak_batch_effect_res$neighbours_original) - 1),
  info = paste(
    "bbknn - biological signal is not regressed out",
    "(weak batch effect)"
  )
)

#### medium batch effects ------------------------------------------------------

medium_batch_effect_res <- assess_bbknn_impact(
  object = sc_object.medium_batch_effect
)

expect_true(
  current = medium_batch_effect_res$kbet_original >
    medium_batch_effect_res$kbet_correct,
  info = paste("bbknn - medium batch effects get regressed out")
)

expect_true(
  current = mean(medium_batch_effect_res$neighbours_corrected) >
    (mean(medium_batch_effect_res$neighbours_original) - 1),
  info = paste(
    "bbknn - biological signal is not regressed out",
    "(medium batch effect)"
  )
)

#### strong batch effects ------------------------------------------------------

strong_batch_effect_res <- assess_bbknn_impact(
  object = sc_object.strong_batch_effect
)

expect_true(
  current = strong_batch_effect_res$kbet_original >
    strong_batch_effect_res$kbet_correct,
  info = paste("bbknn - medium batch effects get regressed out")
)

expect_true(
  current = mean(strong_batch_effect_res$neighbours_corrected) >
    (mean(strong_batch_effect_res$neighbours_original) - 1),
  info = paste(
    "bbknn - biological signal is not regressed out",
    "(medium batch effect)"
  )
)

sc_object.medium_batch_effect[[]]

object = sc_object.medium_batch_effect
batch_column = "batch_index"
embd_to_use = "pca"
k = 10L
no_embd_to_use = NULL
sc_bbknn_params = params_sc_bbknn()
.verbose = TRUE

class(object)

checkmate::assertTRUE(S7::S7_inherits(object, single_cell_exp))
checkmate::qassert(batch_column, "S1")
checkmate::assertChoice(embd_to_use, c("pca"))
checkmate::qassert(no_embd_to_use, c("I1", "0"))
assertScBbknn(sc_bbknn_params)

embd <- switch(embd_to_use, pca = get_pca_factors(object))
# early return
if (is.null(embd)) {
  warning(
    paste(
      "The desired embedding was not found. Please check the parameters.",
      "Returning NULL."
    )
  )

  return(NULL)
}
if (!is.null(no_embd_to_use)) {
  to_take <- min(c(no_embd_to_use, ncol(embd)))
  embd <- embd[, 1:to_take]
}

batch_index <- unlist(object[[batch_column]])

if (!length(levels(factor(batch_index))) > 1) {
  warning("The batch column only has one batch. Returning NULL")
  return(NULL)
}

no_generated_neighbours <- length(levels(factor(batch_index))) *
  sc_bbknn_params$neighbours_within_batch

if (k > no_generated_neighbours) {
  warning(paste(
    "The number of desired neighbours cannot be generated with these BBKNN",
    "Please adopt neighbours_within_batch accordingly"
  ))
}

message("Running BBKNN algorithm")

bbknn_res <- rs_bbknn(
  embd = embd,
  batch_labels = as.integer(batch_index),
  bbknn_params = sc_bbknn_params,
  seed = 42L,
  verbose = .verbose
)

knn_mat <- if (no_generated_neighbours <= k) {
  matrix(data = bbknn_res$distances$indices, nrow = nrow(embd), byrow = TRUE)
} else {
  rs_bbknn_filtering(
    indptr = bbknn_res$distances$indptr,
    indices = bbknn_res$distances$indices,
    k
  )
}

object <- set_knn(object, knn_mat = knn_mat)

message("Generating graph based on BBKNN connectivities.")

sparse_mat <- Matrix::sparseMatrix(
  i = rep(
    seq_along(bbknn_res$connectivities$indptr[-1]),
    diff(bbknn_res$connectivities$indptr)
  ),
  j = bbknn_res$connectivities$indices + 1,
  x = bbknn_res$connectivities$data,
  dims = c(bbknn_res$connectivities$nrow, bbknn_res$connectivities$ncol),
  index1 = TRUE
)

snn_graph <- igraph::graph_from_adjacency_matrix(
  sparse_mat,
  mode = "max",
  weighted = TRUE
)
