# sc subset processing ---------------------------------------------------------

set.seed(123L)

test_temp_dir <- file.path(tempdir(), "subset_processing")
dir.create(test_temp_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot("Test directory does not exist" = dir.exists(test_temp_dir))

## parameters ------------------------------------------------------------------

min_lib_size <- 300L
min_genes_exp <- 45L
min_cells_exp <- 500L
hvg_to_keep <- 30L
no_pcs <- 15L

## set up parent object --------------------------------------------------------

single_cell_test_data <- generate_single_cell_test_data()

sc_qc_param <- params_sc_min_quality(
  min_unique_genes = min_genes_exp,
  min_lib_size = min_lib_size,
  min_cells = min_cells_exp,
  target_size = 1000
)

sc_object <- SingleCells(dir_data = test_temp_dir)

sc_object <- load_r_data(
  object = sc_object,
  counts = single_cell_test_data$counts,
  obs = single_cell_test_data$obs,
  var = single_cell_test_data$var,
  sc_qc_param = sc_qc_param,
  streaming = 0L,
  .verbose = FALSE
)

# tests ------------------------------------------------------------------------

## subset construction ---------------------------------------------------------

target_group <- "cell_type_1"

subset_obj <- SingleCellsSubset(
  sc_object = sc_object,
  grouping_column = "cell_grp",
  group = target_group
)

### shape ----------------------------------------------------------------------

obs_full_filt <- get_sc_obs(sc_object, filtered = TRUE)
n_expected <- sum(obs_full_filt$cell_grp == target_group)

expect_true(
  current = n_expected > 50L,
  info = "sensible number of cells in target group post-QC"
)

expect_equal(
  current = dim(subset_obj),
  target = c(n_expected, nrow(get_sc_var(sc_object))),
  info = "subset dims = (n group cells passing QC, n genes)"
)

expect_equal(
  current = nrow(get_sc_obs(subset_obj)),
  target = n_expected,
  info = "obs restricted to subset cells"
)

expect_true(
  current = all(unlist(subset_obj[["cell_grp"]]) == target_group),
  info = "all obs rows belong to target group"
)

expect_equal(
  current = nrow(get_sc_var(subset_obj)),
  target = nrow(get_sc_var(sc_object)),
  info = "var table unchanged in subset"
)

### subset_to_original ---------------------------------------------------------

s_to_o <- subset_obj@subset_to_original

expect_equal(
  current = length(s_to_o),
  target = n_expected,
  info = "subset_to_original has subset length"
)

expect_equal(
  current = s_to_o,
  target = as.integer(get_sc_obs(subset_obj)$cell_idx),
  info = "subset_to_original matches obs_table$cell_idx"
)

expect_true(
  current = all(s_to_o %in% obs_full_filt$cell_idx),
  info = "subset_to_original entries are valid parent positions"
)

### sc_map restriction ---------------------------------------------------------

sc_map_sub <- get_sc_map(subset_obj)

expect_equal(
  current = sc_map_sub$cell_mapping,
  target = get_sc_map(sc_object)$cell_mapping,
  info = "sc_map$cell_mapping unchanged from parent (preserves rowname lookup)"
)

expect_true(
  current = all((sc_map_sub$cells_to_keep_idx + 1L) %in% s_to_o),
  info = "sc_map$cells_to_keep_idx restricted to subset (0-indexed)"
)

expect_true(
  current = is.null(sc_map_sub$hvg_gene_indices),
  info = "HVG cleared by constructor"
)

### cells_to_keep --------------------------------------------------------------

c_to_k <- get_cells_to_keep(subset_obj)

expect_true(
  current = all((c_to_k + 1L) %in% s_to_o),
  info = "get_cells_to_keep returns 0-indexed original positions within subset"
)

### bad inputs -----------------------------------------------------------------

expect_error(
  current = SingleCellsSubset(sc_object, "cell_grp", "nonexistent"),
  info = "error on empty group"
)

expect_error(
  current = SingleCellsSubset(sc_object, "not_a_column", target_group),
  info = "error on missing grouping column"
)

## indexing translation --------------------------------------------------------

### subset[1:n, ] gives the right cells in the right order ---------------------

n_take <- 5L
counts_via_subset <- as.matrix(subset_obj[seq_len(n_take), ])
target_cell_ids <- get_sc_obs(subset_obj)$cell_id[seq_len(n_take)]
counts_via_parent <- as.matrix(sc_object[target_cell_ids, ])

expect_equivalent(
  current = counts_via_subset,
  target = counts_via_parent,
  info = "subset[1:n, ] returns same counts as parent[matching_cell_ids, ]"
)

expect_equal(
  current = rownames(counts_via_subset),
  target = target_cell_ids,
  info = "subset[1:n, ] rownames are subset cell_ids"
)

### subset[cell_ids, ] resolves correctly -------------------------------------

some_ids <- get_sc_obs(subset_obj)$cell_id[c(2L, 7L, 11L)]

expect_equivalent(
  current = as.matrix(subset_obj[some_ids, ]),
  target = as.matrix(sc_object[some_ids, ]),
  info = "subset[cell_ids, ] matches parent[same_cell_ids, ]"
)

### get_cell_indices contract --------------------------------------------------

test_ids <- get_sc_obs(subset_obj)$cell_id[c(3L, 5L)]

expect_equal(
  current = get_cell_indices(
    subset_obj,
    cell_ids = test_ids,
    rust_index = FALSE
  ),
  target = c(3L, 5L),
  info = "get_cell_indices(rust_index = FALSE) -> 1-indexed subset positions"
)

expect_equal(
  current = get_cell_indices(
    subset_obj,
    cell_ids = test_ids,
    rust_index = TRUE
  ),
  target = as.integer(s_to_o[c(3L, 5L)] - 1L),
  info = "get_cell_indices(rust_index = TRUE) -> 0-indexed original positions"
)

### get_gene_indices is unchanged (no gene subsetting) ------------------------

gene_test <- get_sc_var(subset_obj)$gene_id[c(4L, 10L)]

expect_equal(
  current = get_gene_indices(
    subset_obj,
    gene_ids = gene_test,
    rust_index = FALSE
  ),
  target = c(4L, 10L),
  info = "get_gene_indices unchanged in subset (gene space identical)"
)

### bounds checking ------------------------------------------------------------

expect_error(
  current = get_sc_counts(subset_obj, cell_indices = n_expected + 1L),
  info = "out-of-bounds subset position errors"
)

## obs setters -----------------------------------------------------------------

new_label <- rep("foo", nrow(get_sc_obs(subset_obj)))
subset_obj[["my_label"]] <- new_label

expect_true(
  current = all(unlist(subset_obj[["my_label"]]) == "foo"),
  info = "[[<- adds an obs column on subset"
)

## HVG on subset ---------------------------------------------------------------

subset_obj <- find_hvg_sc(subset_obj, hvg_no = hvg_to_keep, .verbose = FALSE)

hvg_subset <- get_hvg(subset_obj)

expect_true(
  current = checkmate::qtest(hvg_subset, sprintf("I%i[0,)", hvg_to_keep)),
  info = "subset HVG: correct length, 0-indexed"
)

expect_equal(
  current = get_sc_map(subset_obj)$hvg_gene_indices,
  target = hvg_subset,
  info = "HVG stored in sc_map (0-indexed) and round-trips via get_hvg"
)

### subset HVG diverges from full HVG ------------------------------------------

# Full data HVG should be dominated by between-cell-type marker variance
# (gene indices 0-29 0-indexed). Subset HVG should NOT preferentially include
# the target cell type's own markers (0-9 for cell_type_1), since they are
# constitutive within the subset and have low variance.
sc_for_hvg <- find_hvg_sc(sc_object, hvg_no = hvg_to_keep, .verbose = FALSE)
hvg_full <- get_hvg(sc_for_hvg)

ct1_markers_0idx <- 0:9L

expect_true(
  current = length(intersect(hvg_subset, ct1_markers_0idx)) <
    length(intersect(hvg_full, ct1_markers_0idx)),
  info = "ct1 markers under-represented in subset HVG vs full HVG"
)

expect_false(
  current = identical(sort(hvg_subset), sort(hvg_full)),
  info = "subset HVG differs from full HVG"
)

## PCA on subset ---------------------------------------------------------------

subset_obj <- calculate_pca_sc(
  object = subset_obj,
  no_pcs = no_pcs,
  pca_params = params_sc_pca(),
  .verbose = FALSE
)

pca_subset <- get_pca_factors(subset_obj)

expect_true(
  current = checkmate::testMatrix(
    x = pca_subset,
    mode = "numeric",
    nrows = n_expected,
    ncols = no_pcs,
    row.names = "named",
    col.names = "named"
  ),
  info = "subset PCA factors: shape (n_subset, no_pcs), named"
)

expect_equal(
  current = rownames(pca_subset),
  target = get_sc_obs(subset_obj)$cell_id,
  info = "PCA factors rownames match subset obs cell_ids in order"
)

### PCA loadings rownames use subset HVG --------------------------------------

pca_loadings <- get_pca_loadings(subset_obj)

expect_equal(
  current = rownames(pca_loadings),
  target = get_sc_var(subset_obj)$gene_id[hvg_subset + 1L],
  info = "PCA loadings rownames = subset HVG gene_ids"
)

### subset PCA matches direct Rust call with same cells + HVG -----------------

zeallot::`%<-%`(
  c(pca_direct, ., ., .),
  rs_sc_pca(
    f_path_gene = bixverse:::get_rust_count_gene_f_path(sc_object),
    f_path_cell = bixverse:::get_rust_count_cell_f_path(sc_object),
    no_pcs = no_pcs,
    pca_params = params_sc_pca(),
    cell_indices = get_cells_to_keep(subset_obj),
    gene_indices = hvg_subset,
    seed = 42L,
    return_scaled = FALSE,
    verbose = 0L
  )
)

expect_equal(
  current = abs(diag(cor(pca_subset, pca_direct))),
  target = rep(1, no_pcs),
  tolerance = 1e-8,
  info = "subset PCA == direct Rust call with same cells + same HVG"
)

## fast clustering on subset ---------------------------------------------------

fc_res <- fast_cluster_sc(
  object = subset_obj,
  .verbose = FALSE
)

expect_true(
  current = inherits(fc_res, "SingleCellFastClusters"),
  info = "fast_cluster_sc returns SingleCellFastClusters"
)

expect_equal(
  current = nrow(fc_res$memberships),
  target = n_expected,
  info = "memberships rows match subset cell count"
)

expect_true(
  current = all(fc_res$memberships$cell_idx %in% s_to_o),
  info = "memberships$cell_idx in original index space (matches obs_table$cell_idx)"
)

## neighbours + clustering on subset (via ScOrMc) ------------------------------

subset_obj <- find_neighbours_sc(subset_obj, .verbose = FALSE)

expect_equal(
  current = nrow(get_knn_mat(subset_obj)),
  target = n_expected,
  info = "kNN matrix row count matches subset cells"
)

expect_true(
  current = inherits(get_snn_graph(subset_obj), "igraph"),
  info = "sNN graph stored on subset"
)

subset_obj <- find_clusters_sc(subset_obj, name = "leiden_clustering")

expect_true(
  current = "leiden_clustering" %in% colnames(get_sc_obs(subset_obj)),
  info = "find_clusters_sc writes column to subset obs"
)

expect_equal(
  current = length(unlist(subset_obj[["leiden_clustering"]])),
  target = n_expected,
  info = "leiden column has one entry per subset cell"
)

## state isolation: parent untouched by subset ops -----------------------------

# subset operations should not leak back into the parent.
expect_true(
  current = is.null(get_sc_map(sc_object)$hvg_gene_indices) ||
    !identical(get_sc_map(sc_object)$hvg_gene_indices, hvg_subset),
  info = "parent HVG not mutated by subset find_hvg_sc"
)

expect_true(
  current = is.null(get_pca_factors(sc_object)),
  info = "parent PCA cache not mutated by subset calculate_pca_sc"
)

# cleanup ----------------------------------------------------------------------

on.exit(unlink(test_temp_dir, recursive = TRUE, force = TRUE), add = TRUE)
