# sc subset batch correction ---------------------------------------------------

library(magrittr)

test_temp_dir <- file.path(tempdir(), "sc_subset_batch_corr")
dir.create(test_temp_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot("Test directory does not exist" = dir.exists(test_temp_dir))

## testing parameters ----------------------------------------------------------

min_lib_size <- 300L
min_genes_exp <- 45L
min_cells_exp <- 450L
hvg_to_keep <- 30L
no_pcs <- 15L

cell_markers <- list(
  cell_type_1 = list(marker_genes = 0:8L),
  cell_type_2 = list(marker_genes = 9:19L),
  cell_type_3 = list(marker_genes = 20:29L),
  cell_type_4 = list(marker_genes = 30:44L)
)

target_group <- "cell_type_1"

## helpers ---------------------------------------------------------------------

# Build a SingleCellsSubset with HVG + PCA + kNN ready for batch testing.
build_subset_for_batch <- function(strength, dir_data) {
  test_data <- generate_single_cell_test_data(
    syn_data_params = params_sc_synthetic_data(
      n_cells = 900L,
      marker_genes = cell_markers,
      n_batches = 3L,
      batch_effect_strength = strength
    )
  )

  parent <- SingleCells(dir_data = dir_data)
  parent <- load_r_data(
    object = parent,
    counts = test_data$counts,
    obs = test_data$obs,
    var = test_data$var,
    sc_qc_param = params_sc_min_quality(
      min_unique_genes = min_genes_exp,
      min_lib_size = min_lib_size,
      min_cells = min_cells_exp
    ),
    streaming = 0L,
    .verbose = FALSE
  )

  subset_obj <- SingleCellsSubset(
    sc_object = parent,
    grouping_column = "cell_grp",
    group = target_group
  )

  subset_obj <- find_hvg_sc(subset_obj, hvg_no = hvg_to_keep, .verbose = FALSE)
  subset_obj <- calculate_pca_sc(subset_obj, no_pcs = no_pcs, .verbose = FALSE)
  subset_obj <- find_neighbours_sc(
    subset_obj,
    neighbours_params = params_sc_neighbours(knn = list(k = 10L)),
    .verbose = FALSE
  )

  subset_obj
}

# kBET rejection rate on the current kNN of a subset.
kbet_rate <- function(subset_obj) {
  batch <- as.integer(unlist(subset_obj[["batch_index"]]))
  res <- rs_kbet(
    knn_mat = get_knn_mat(subset_obj),
    batch_vector = batch,
    verbose = FALSE
  )
  sum(res$pval <= 0.05) / length(res$pval)
}

## pre-processing --------------------------------------------------------------

dir_data <- list("weak" = NULL, "medium" = NULL, "strong" = NULL)
dir_data <- purrr::imap(dir_data, \(elem, name) {
  final_path <- file.path(test_temp_dir, name)
  dir.create(final_path, showWarnings = FALSE)
  final_path
})

subset.weak <- build_subset_for_batch("weak", dir_data$weak)
subset.medium <- build_subset_for_batch("medium", dir_data$medium)
subset.strong <- build_subset_for_batch("strong", dir_data$strong)

## sanity ----------------------------------------------------------------------

expect_true(
  current = nrow(get_sc_obs(subset.medium)) > 100L,
  info = "subset has sensible number of cells post-QC"
)

expect_true(
  current = length(unique(unlist(subset.medium[["batch_index"]]))) == 3L,
  info = "all three batches represented in the subset"
)

expect_true(
  current = all(unlist(subset.medium[["cell_grp"]]) == target_group),
  info = "subset obs restricted to target group"
)

# tests ------------------------------------------------------------------------

## metrics ---------------------------------------------------------------------

### kbet scores ----------------------------------------------------------------

kbet.weak <- calculate_kbet_sc(
  object = subset.weak,
  batch_column = "batch_index",
  .verbose = FALSE
)

kbet.medium <- calculate_kbet_sc(
  object = subset.medium,
  batch_column = "batch_index",
  .verbose = FALSE
)

kbet.strong <- calculate_kbet_sc(
  object = subset.strong,
  batch_column = "batch_index",
  .verbose = FALSE
)

expect_true(
  current = is.logical(kbet.weak$significant_tests),
  info = "subset kbet - returns booleans where expected"
)

expect_true(
  current = checkmate::qtest(kbet.weak$p_values, "N[0, 1]"),
  info = "subset kbet - p-values are expected type and range"
)

expect_true(
  current = kbet.weak$kbet_score < kbet.medium$kbet_score,
  info = "subset kbet - weak < medium"
)

expect_true(
  current = kbet.medium$kbet_score < kbet.strong$kbet_score,
  info = "subset kbet - medium < strong"
)

### silhouette scores ----------------------------------------------------------

asw.weak <- calculate_batch_asw_sc(
  object = subset.weak,
  batch_column = "batch_index",
  .verbose = FALSE
)

asw.medium <- calculate_batch_asw_sc(
  object = subset.medium,
  batch_column = "batch_index",
  .verbose = FALSE
)

asw.strong <- calculate_batch_asw_sc(
  object = subset.strong,
  batch_column = "batch_index",
  .verbose = FALSE
)

expect_true(
  current = checkmate::qtest(asw.weak$per_cell, "N+[-1, 1]"),
  info = "subset asw - per_cell in [-1, 1]"
)

expect_true(
  current = checkmate::qtest(asw.weak$mean_asw, "N1[-1, 1]"),
  info = "subset asw - mean in range"
)

expect_true(
  current = asw.weak$mean_asw < asw.medium$mean_asw,
  info = "subset asw - weak < medium"
)

expect_true(
  current = asw.medium$mean_asw < asw.strong$mean_asw,
  info = "subset asw - medium < strong"
)

### lisi scores ----------------------------------------------------------------

lisi.weak <- calculate_batch_lisi_sc(
  object = subset.weak,
  batch_column = "batch_index",
  .verbose = FALSE
)

lisi.medium <- calculate_batch_lisi_sc(
  object = subset.medium,
  batch_column = "batch_index",
  .verbose = FALSE
)

lisi.strong <- calculate_batch_lisi_sc(
  object = subset.strong,
  batch_column = "batch_index",
  .verbose = FALSE
)

expect_true(
  current = checkmate::qtest(lisi.weak$per_cell, "N+"),
  info = "subset lisi - per_cell numeric"
)

# higher LISI = better mixing
expect_true(
  current = lisi.medium$mean_lisi < lisi.weak$mean_lisi,
  info = "subset lisi - medium < weak"
)

expect_true(
  current = lisi.strong$mean_lisi < lisi.medium$mean_lisi,
  info = "subset lisi - strong < medium"
)

## hvg batch aware -------------------------------------------------------------

hvg_union <- find_hvg_batch_aware_sc(
  subset.strong,
  hvg_no = 30L,
  batch_column = "batch_index",
  gene_comb_method = "union",
  .verbose = FALSE
)

hvg_avg <- find_hvg_batch_aware_sc(
  subset.strong,
  hvg_no = 30L,
  batch_column = "batch_index",
  gene_comb_method = "average",
  .verbose = FALSE
)

hvg_intersect <- find_hvg_batch_aware_sc(
  subset.strong,
  hvg_no = 30L,
  batch_column = "batch_index",
  gene_comb_method = "intersection",
  .verbose = FALSE
)

expect_true(
  current = checkmate::qtest(hvg_union$hvg_gene_idx, "I+"),
  info = "subset hvg batch aware - indices (union)"
)

expect_true(
  current = checkmate::qtest(hvg_union$hvg_genes, "S+"),
  info = "subset hvg batch aware - gene names (union)"
)

expect_true(
  current = checkmate::testDataTable(hvg_union$hvg_data),
  info = "subset hvg batch aware - data.table (union)"
)

expect_true(
  current = length(hvg_union$hvg_genes) > length(hvg_avg$hvg_genes),
  info = "subset hvg batch aware - union > average"
)

expect_true(
  current = length(hvg_avg$hvg_genes) > length(hvg_intersect$hvg_genes),
  info = "subset hvg batch aware - average > intersection"
)

expect_true(
  current = length(hvg_avg$hvg_genes) == 30L,
  info = "subset hvg batch aware - average returns hvg_no"
)

## bbknn -----------------------------------------------------------------------

# warning when overwriting existing kNN
expect_warning(
  current = bbknn_sc(
    object = subset.medium,
    batch_column = "batch_index",
    .verbose = FALSE
  ),
  info = "subset bbknn - warning when overwriting existing kNN"
)

# warning when requested neighbours > generated
expect_warning(
  current = bbknn_sc(
    object = remove_knn(subset.medium),
    batch_column = "batch_index",
    no_neighbours_to_keep = 25L,
    .verbose = FALSE
  ),
  info = "subset bbknn - warning when too many neighbours requested"
)

kbet_pre_bbknn <- calculate_kbet_sc(
  object = subset.medium,
  batch_column = "batch_index",
  .verbose = FALSE
)$kbet_score

subset_bbknn <- bbknn_sc(
  object = remove_knn(subset.medium),
  batch_column = "batch_index",
  no_neighbours_to_keep = 9L,
  .verbose = FALSE
)

kbet_post_bbknn <- calculate_kbet_sc(
  object = subset_bbknn,
  batch_column = "batch_index",
  .verbose = FALSE
)$kbet_score

expect_true(
  current = kbet_pre_bbknn > kbet_post_bbknn,
  info = "subset bbknn - batch effect reduced"
)

expect_true(
  current = checkmate::testMatrix(
    get_knn_mat(subset_bbknn),
    mode = "integer",
    nrow = nrow(get_sc_obs(subset.medium)),
    ncol = 9L
  ),
  info = "subset bbknn - kNN matrix correct type and dimensions"
)

expect_true(
  current = checkmate::testClass(get_snn_graph(subset_bbknn), "igraph"),
  info = "subset bbknn - sNN graph generated"
)

## fastmnn ---------------------------------------------------------------------

# Assess kBET before/after fastMNN on a subset.
assess_subset_fastmnn <- function(subset_obj) {
  hvg_ba <- find_hvg_batch_aware_sc(
    subset_obj,
    hvg_no = 30L,
    batch_column = "batch_index",
    gene_comb_method = "union",
    .verbose = FALSE
  )
  before <- kbet_rate(subset_obj)
  corrected <- fast_mnn_sc(
    object = subset_obj,
    batch_column = "batch_index",
    batch_hvg_genes = hvg_ba$hvg_gene_idx,
    fastmnn_params = params_sc_fastmnn(no_pcs = 10L, knn = list(k = 5L)),
    .verbose = FALSE
  )
  corrected <- find_neighbours_sc(
    corrected,
    embd_to_use = "mnn",
    neighbours_params = params_sc_neighbours(knn = list(k = 10L)),
    .verbose = FALSE
  )
  list(before = before, after = kbet_rate(corrected), corrected = corrected)
}

### weak batch effects ---------------------------------------------------------

mnn.weak <- assess_subset_fastmnn(subset.weak)

expect_true(
  current = mnn.weak$before > mnn.weak$after,
  info = "subset fastmnn - weak batch effect reduced"
)

### medium batch effects -------------------------------------------------------

mnn.medium <- assess_subset_fastmnn(subset.medium)

expect_true(
  current = mnn.medium$before > mnn.medium$after,
  info = "subset fastmnn - medium batch effect reduced"
)

### strong batch effects -------------------------------------------------------

mnn.strong <- assess_subset_fastmnn(subset.strong)

expect_true(
  current = mnn.strong$before > mnn.strong$after,
  info = "subset fastmnn - strong batch effect reduced"
)

### embedding shape ------------------------------------------------------------

expect_true(
  current = checkmate::testMatrix(
    get_embedding(mnn.medium$corrected, "mnn"),
    mode = "numeric",
    nrows = nrow(get_sc_obs(subset.medium))
  ),
  info = "subset fastmnn - mnn embedding has correct row count"
)

## harmony ---------------------------------------------------------------------

# Assess kBET before/after Harmony (v1 or v2) on a subset.
assess_subset_harmony <- function(
  subset_obj,
  harmony_fn,
  params_fn,
  embd_name
) {
  before <- kbet_rate(subset_obj)
  corrected <- harmony_fn(
    object = subset_obj,
    batch_column = "batch_index",
    harmony_params = params_fn(),
    .verbose = FALSE
  )
  corrected <- find_neighbours_sc(
    corrected,
    embd_to_use = embd_name,
    neighbours_params = params_sc_neighbours(knn = list(k = 10L)),
    .verbose = FALSE
  )
  list(before = before, after = kbet_rate(corrected), corrected = corrected)
}

### weak batch effects ---------------------------------------------------------

#### harmony v1 ----------------------------------------------------------------

h1.weak <- assess_subset_harmony(
  subset.weak,
  harmony_sc,
  params_sc_harmony,
  "harmony"
)

expect_true(
  current = h1.weak$before > h1.weak$after,
  info = "subset harmony - weak batch effect reduced"
)

#### harmony v2 ----------------------------------------------------------------

h2.weak <- assess_subset_harmony(
  subset.weak,
  harmony_v2_sc,
  params_sc_harmony_v2,
  "harmony_v2"
)

expect_true(
  current = h2.weak$before > h2.weak$after,
  info = "subset harmony_v2 - weak batch effect reduced"
)

### medium batch effects -------------------------------------------------------

#### harmony v1 ----------------------------------------------------------------

h1.medium <- assess_subset_harmony(
  subset.medium,
  harmony_sc,
  params_sc_harmony,
  "harmony"
)

expect_true(
  current = h1.medium$before > h1.medium$after,
  info = "subset harmony - medium batch effect reduced"
)

#### harmony v2 ----------------------------------------------------------------

h2.medium <- assess_subset_harmony(
  subset.medium,
  harmony_v2_sc,
  params_sc_harmony_v2,
  "harmony_v2"
)

expect_true(
  current = h2.medium$before > h2.medium$after,
  info = "subset harmony_v2 - medium batch effect reduced"
)

### strong batch effects -------------------------------------------------------

#### harmony v1 ----------------------------------------------------------------

h1.strong <- assess_subset_harmony(
  subset.strong,
  harmony_sc,
  params_sc_harmony,
  "harmony"
)

expect_true(
  current = h1.strong$before > h1.strong$after,
  info = "subset harmony - strong batch effect reduced"
)

#### harmony v2 ----------------------------------------------------------------

h2.strong <- assess_subset_harmony(
  subset.strong,
  harmony_v2_sc,
  params_sc_harmony_v2,
  "harmony_v2"
)

expect_true(
  current = h2.strong$before > h2.strong$after,
  info = "subset harmony_v2 - strong batch effect reduced"
)

### embedding shape ------------------------------------------------------------

expect_true(
  current = checkmate::testMatrix(
    get_embedding(h1.medium$corrected, "harmony"),
    mode = "numeric",
    nrows = nrow(get_sc_obs(subset.medium))
  ),
  info = "subset harmony - embedding has correct row count"
)

expect_true(
  current = checkmate::testMatrix(
    get_embedding(h2.medium$corrected, "harmony_v2"),
    mode = "numeric",
    nrows = nrow(get_sc_obs(subset.medium))
  ),
  info = "subset harmony_v2 - embedding has correct row count"
)

# clean up ---------------------------------------------------------------------

on.exit(unlink(test_temp_dir, recursive = TRUE, force = TRUE), add = TRUE)
