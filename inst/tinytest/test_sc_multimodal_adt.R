# sc multi-modal (adt) ----------------------------------------------------------

set.seed(123L)

test_temp_dir_mm <- file.path(tempdir(), "processing_mm")
dir.create(test_temp_dir_mm, recursive = TRUE, showWarnings = FALSE)
stopifnot("Test directory does not exist" = dir.exists(test_temp_dir_mm))

## testing parameters ----------------------------------------------------------

min_lib_size <- 300L
min_genes_exp <- 45L
min_cells_exp <- 500L
hvg_to_keep <- 30L
no_pcs <- 15L
no_pcs_adt <- 10L

## synthetic test data ---------------------------------------------------------

rna <- generate_single_cell_test_data()
adt <- generate_single_cell_test_data_adt()

# the two modalities must describe the same barcodes for multi-modal to align
expect_true(
  current = all(rna$obs$cell_id %in% rownames(adt$counts)),
  info = "rna and adt synthetic data share barcodes"
)

## underlying class + rna ingestion --------------------------------------------

sc_mm <- SingleCellsMultiModal(dir_data = test_temp_dir_mm)

sc_mm <- load_r_data(
  object = sc_mm,
  counts = rna$counts,
  obs = rna$obs,
  var = rna$var,
  sc_qc_param = params_sc_min_quality(
    min_unique_genes = min_genes_exp,
    min_lib_size = min_lib_size,
    min_cells = min_cells_exp
  ),
  streaming = 0L,
  .verbose = FALSE
)

expect_true(
  current = checkmate::qtest(get_cells_to_keep(sc_mm), "I+"),
  info = "mm: cells_to_keep set after rna ingestion"
)

## add adt counts --------------------------------------------------------------

# kept set at the time the ADTCounts is frozen
keep_full <- get_cell_names(sc_mm, filtered = TRUE)
adt_raw_kept <- adt$counts[keep_full, ]

sc_mm <- add_adt_counts_sc(
  sc_mm,
  adt_counts = adt$counts,
  method = "clr"
)

### obs / var population --------------------------------------------------------

obs_after <- get_sc_obs(sc_mm, filtered = TRUE)

expect_true(
  current = all(c("adt_nnz", "adt_lib_size") %in% names(obs_after)),
  info = "mm: adt sample info joined into obs"
)

# adt_lib_size must equal the raw rowSums of the kept cells (exact)
lib_ref <- Matrix::rowSums(adt_raw_kept)
lib_obs <- obs_after$adt_lib_size[match(names(lib_ref), obs_after$cell_id)]
expect_equivalent(
  current = lib_obs,
  target = as.numeric(lib_ref),
  tolerance = 1e-8,
  info = "mm: adt_lib_size matches manual rowSums"
)

adt_var <- get_sc_var(sc_mm, modality = "adt")

expect_true(
  current = nrow(adt_var) == adt$counts %>% ncol() &&
    "feature_id" %in% names(adt_var),
  info = "mm: var_adt populated with one row per protein"
)

expect_true(
  current = all(get_adt_names(sc_mm) == colnames(adt$counts)),
  info = "mm: adt feature names retrievable"
)

## normalisation ---------------------------------------------------------------

### biological signal -----------------------------------------------------------

# stored norm uses the default clean_clr_counts = TRUE; markers must still be
# elevated in their own cell type.
adt_norm <- get_sc_counts(sc_mm, modality = "adt", assay = "norm")
grp <- obs_after$cell_grp[match(rownames(adt_norm), obs_after$cell_id)]

m1 <- tapply(adt_norm[, "protein_01"], grp, mean)
m2 <- tapply(adt_norm[, "protein_04"], grp, mean)
m3 <- tapply(adt_norm[, "protein_07"], grp, mean)

expect_true(
  current = m1["cell_type_1"] > m1["cell_type_2"] &&
    m1["cell_type_1"] > m1["cell_type_3"],
  info = "mm: marker protein elevated in its own cell type (1)"
)

expect_true(
  current = m2["cell_type_2"] > m2["cell_type_1"] &&
    m2["cell_type_2"] > m2["cell_type_3"],
  info = "mm: marker protein elevated in its own cell type (2)"
)

expect_true(
  current = m3["cell_type_3"] > m3["cell_type_1"] &&
    m3["cell_type_3"] > m3["cell_type_2"],
  info = "mm: marker protein elevated in its own cell type (3)"
)

## filtering reduces returned counts -------------------------------------------

n_full <- length(get_cells_to_keep(sc_mm))
expect_true(
  current = nrow(get_sc_counts(sc_mm, modality = "adt", assay = "norm")) ==
    n_full,
  info = "mm: adt accessor returns the full kept set"
)

keep_subset <- keep_full[1:400]
sc_mm <- set_cells_to_keep(sc_mm, keep_subset)

adt_sub <- get_sc_counts(sc_mm, modality = "adt", assay = "norm")
expect_true(
  current = nrow(adt_sub) == length(keep_subset),
  info = "mm: adt accessor shrinks with a tighter cells_to_keep"
)
expect_true(
  current = all(rownames(adt_sub) == keep_subset),
  info = "mm: adt accessor returns the current kept barcodes in order"
)

# restore the fuller kept set for downstream analysis (also tests name-keyed
# re-expansion against the frozen matrix)
sc_mm <- set_cells_to_keep(sc_mm, keep_full)
expect_true(
  current = nrow(get_sc_counts(sc_mm, modality = "adt", assay = "norm")) ==
    n_full,
  info = "mm: adt accessor re-expands when filter is relaxed"
)

## pca on adt ------------------------------------------------------------------

adt_norm <- get_sc_counts(sc_mm, modality = "adt", assay = "norm")
pca_r_adt <- prcomp(as.matrix(adt_norm), scale. = TRUE)

sc_mm <- calculate_pca_adt_sc(
  object = sc_mm,
  no_pcs = no_pcs_adt,
  randomised_svd = FALSE
)

adt_factors <- get_pca_factors(sc_mm, modality = "adt")

expect_true(
  current = checkmate::testMatrix(
    adt_factors,
    mode = "numeric",
    row.names = "named",
    col.names = "named"
  ),
  info = "mm: adt pca factors are a named matrix"
)

expect_true(
  current = all.equal(
    abs(diag(cor(
      adt_factors[, 1:no_pcs_adt],
      pca_r_adt$x[, 1:no_pcs_adt]
    ))),
    rep(1, no_pcs_adt),
    tolerance = 1e-6
  ),
  info = "mm: adt pca matches prcomp"
)

## knn on adt ------------------------------------------------------------------

sc_mm <- find_neighbours_sc(
  sc_mm,
  embd_to_use = "pca",
  modality = "adt",
  .verbose = FALSE
)

adt_knn <- get_knn_mat(sc_mm, modality = "adt")
expect_true(
  current = nrow(adt_knn) == length(get_cells_to_keep(sc_mm)),
  info = "mm: adt knn has one row per kept cell"
)
expect_true(
  current = checkmate::testClass(
    get_snn_graph(sc_mm, modality = "adt"),
    "igraph"
  ),
  info = "mm: adt snn graph returned as igraph"
)

## wnn -------------------------------------------------------------------------

# wnn needs a pca embedding on both modalities
sc_mm <- find_hvg_sc(sc_mm, hvg_no = hvg_to_keep, .verbose = FALSE)
sc_mm <- calculate_pca_sc(
  sc_mm,
  no_pcs = no_pcs,
  pca_params = params_sc_pca(),
  .verbose = FALSE
)

sc_mm <- generate_wnn_graph_sc(
  sc_mm,
  modality_1 = "rna",
  modality_2 = "adt",
  .verbose = FALSE
)

expect_true(
  current = !is.null(get_knn_obj(sc_mm, modality = "wnn")),
  info = "mm: wnn knn stored and retrievable"
)
expect_true(
  current = checkmate::testClass(
    get_snn_graph(sc_mm, modality = "wnn"),
    "igraph"
  ),
  info = "mm: wnn snn graph returned as igraph"
)

# the integration test: clustering on the fused graph recovers cell types
sc_mm <- find_clusters_sc(
  sc_mm,
  modality = "wnn",
  name = "wnn_clustering"
)

cell_grps <- unlist(sc_mm[["cell_grp"]])
wnn_clusters <- unlist(sc_mm[["wnn_clustering"]])

f1_scores <- f1_score_confusion_mat(cell_grps, wnn_clusters)
expect_true(
  current = all(f1_scores > 0.95),
  info = "mm: wnn clustering recovers the cell groups"
)

# clean up ---------------------------------------------------------------------

on.exit(unlink(test_temp_dir_mm, recursive = TRUE, force = TRUE), add = TRUE)
