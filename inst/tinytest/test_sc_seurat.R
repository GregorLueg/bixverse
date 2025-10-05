# comparisons with seurat ------------------------------------------------------

if (!requireNamespace("Seurat", quietly = TRUE)) {
  exit_file("Seurat not available")
}

## synthetic data --------------------------------------------------------------

# thresholds
# absurd numbers, but this is due to the synthetic data
min_lib_size <- 300L
min_genes_exp <- 45L
min_cells_exp <- 500L

single_cell_test_data <- generate_single_cell_test_data()

## seurat tests ----------------------------------------------------------------

counts_transposed <- Matrix::t(single_cell_test_data$counts)

# comparison to seurat... conversion, similar results, etc. pp.
# key differences are likely due to f16 vs f64 in terms of count storage
# for norm counts

### transformation from seurat to bixverse single cell -----------------------

seurat_obj <- suppressWarnings(Seurat::CreateSeuratObject(
  counts = counts_transposed,
  project = "test",
  meta.data = single_cell_test_data$obs,
  min.cells = min_cells_exp
))

sc_object <- suppressWarnings(single_cell_exp(dir_data = tempdir()))

sc_object <- load_seurat(
  object = sc_object,
  seurat = seurat_obj,
  sc_qc_param = params_sc_min_quality(
    min_unique_genes = min_genes_exp,
    min_lib_size = min_lib_size,
    min_cells = min_cells_exp
  ),
  streaming = FALSE,
  .verbose = FALSE
)

# subset here... Seurat is weird as f--k
# also, otherwise dual filtering is happening
seurat_obj <- subset(
  seurat_obj,
  subset = nCount_RNA >= min_lib_size & nFeature_RNA >= min_genes_exp
)

expect_equal(
  current = get_cell_names(sc_object),
  target = colnames(seurat_obj),
  info = "seurat <> bixverse - same cell names"
)

expect_equal(
  current = get_gene_names(sc_object),
  target = rownames(seurat_obj),
  info = "seurat <> bixverse - same gene names"
)

expect_equivalent(
  current = unlist(sc_object[["lib_size"]]),
  target = unlist(seurat_obj$nCount_RNA),
  info = "seurat <> bixverse - same lib size"
)

expect_equivalent(
  current = unlist(sc_object[["nnz"]]),
  target = unlist(seurat_obj$nFeature_RNA),
  info = "seurat <> bixverse - same lib size"
)

### count dimensions ---------------------------------------------------------

seurat_counts <- Seurat::GetAssayData(seurat_obj, layer = "counts")
sc_object_counts <- sc_object[,, return_format = "cell"]

# seurat stores genes x cells in CSC; bixverse stores cells x genes in CSR
# it's the same, but for dimensions and naming of the variables in R
expect_equal(
  current = sc_object_counts@p,
  target = seurat_counts@p,
  info = "seurat <> bixverse - raw counts - same indptr",
)

expect_equal(
  current = sc_object_counts@j,
  target = seurat_counts@i,
  info = "seurat <> bixverse - raw counts - same indices",
)

expect_equal(
  current = sc_object_counts@x,
  target = seurat_counts@x,
  info = "seurat <> bixverse - raw counts - same data",
)

### normalisations -----------------------------------------------------------

seurat_obj <- Seurat::NormalizeData(seurat_obj, verbose = FALSE)

seurat_norm_counts <- Seurat::GetAssayData(seurat_obj, layer = "data")

sc_object_norm_counts <- sc_object[,, assay = "norm", return_format = "cell"]

expect_equal(
  current = seurat_norm_counts@p,
  target = seurat_norm_counts@p,
  info = "seurat <> bixverse - norm counts - same indptr",
)

expect_equal(
  current = sc_object_norm_counts@j,
  target = seurat_norm_counts@i,
  info = "seurat <> bixverse - norm counts - same indices",
)

expect_equal(
  current = sc_object_norm_counts@x,
  target = seurat_norm_counts@x,
  info = "seurat <> bixverse - same counts",
  tolerance = 1e-3 # f16 vs f64!
)

### gene set proportions -----------------------------------------------------

#### calculation -------------------------------------------------------------

seurat_obj[["gs_1"]] <- Seurat::PercentageFeatureSet(
  object = seurat_obj,
  features = c("gene-001", "gene-002", "gene-003", "gene-004")
)
seurat_obj[["gs_2"]] <- Seurat::PercentageFeatureSet(
  object = seurat_obj,
  features = c("gene-096", "gene-097", "gene-100")
)

sc_object <- gene_set_proportions_sc(
  object = sc_object,
  gene_set_list = list(
    gs_1 = c("gene-001", "gene-002", "gene-003", "gene-004"),
    gs_2 = c("gene-096", "gene-097", "gene-100")
  ),
  .verbose = FALSE
)

expect_equivalent(
  current = unlist(sc_object[["gs_1"]]) * 100,
  target = unlist(seurat_obj[["gs_1"]]),
  info = "seurat <> bixverse - same gene proportion calculation",
  tolerance = 1e-7
)

#### filtering ---------------------------------------------------------------

threshold <- 0.05
cells_to_keep <- sc_object[[]][gs_2 < threshold, cell_id]

expect_true(
  current = length(cells_to_keep) > 600,
  info = "sensible cell filtering based on the threshold"
)

sc_object <- set_cell_to_keep(sc_object, cells_to_keep)

seurat_obj <- subset(seurat_obj, subset = gs_2 < threshold * 100)

expect_equal(
  current = get_cell_names(sc_object, filtered = TRUE),
  target = colnames(seurat_obj),
  info = "seurat <> bixverse - same cell names after filtering",
)

### hvg ----------------------------------------------------------------------

hvgs_to_keep <- 30L

seurat_obj <- Seurat::FindVariableFeatures(
  seurat_obj,
  selection.method = "vst",
  nfeatures = hvgs_to_keep,
  verbose = FALSE
)

sc_object <- find_hvg_sc(
  object = sc_object,
  hvg_no = hvgs_to_keep,
  .verbose = FALSE
)

expect_true(
  current = length(intersect(
    Seurat::VariableFeatures(seurat_obj),
    get_gene_names(sc_object)[get_hvg(sc_object) + 1]
  )) ==
    hvgs_to_keep,
  info = "seurat <> bixverse - same hvg identified"
)

### pca ----------------------------------------------------------------------

seurat_obj <- Seurat::ScaleData(
  seurat_obj,
  features = Seurat::VariableFeatures(seurat_obj),
  verbose = FALSE
)

seurat_obj <- Seurat::RunPCA(
  seurat_obj,
  npcs = 10L,
  features = Seurat::VariableFeatures(object = seurat_obj),
  verbose = FALSE,
  weight.by.var = TRUE
)

sc_object <- calculate_pca_sc(
  object = sc_object,
  no_pcs = 10L,
  randomised_svd = FALSE,
  .verbose = FALSE
)

expect_true(
  current = all(
    abs(diag(cor(
      Seurat::Embeddings(seurat_obj),
      get_pca_factors(sc_object)
    ))) >=
      0.99
  ),
  info = "seurat <> bixverse - PCAs highly correlated"
)

sc_object <- calculate_pca_sc(
  object = sc_object,
  no_pcs = 10L,
  randomised_svd = TRUE,
  .verbose = FALSE
)

expect_true(
  current = all(
    abs(diag(cor(
      Seurat::Embeddings(seurat_obj),
      get_pca_factors(sc_object)
    ))) >=
      0.99
  ),
  info = "seurat <> bixverse - PCAs highly correlated (randomised SVD)"
)

### neighbours ---------------------------------------------------------------

# seurat part
seurat_obj <- Seurat::FindNeighbors(
  seurat_obj,
  dims = 1:10,
  k.param = 15L,
  verbose = FALSE
)

seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = 1, verbose = FALSE)

#bixverse part
sc_object <- find_neighbours_sc(
  sc_object,
  neighbours_params = params_sc_neighbours(knn_algorithm = "annoy"),
  .verbose = FALSE
)

sc_object <- find_clusters_sc(sc_object)

seurat_clusters <- seurat_obj$seurat_clusters

bixverse_clusters <- unlist(sc_object[["leiden_clustering"]])

f1_scores <- f1_score_confusion_mat(seurat_clusters, bixverse_clusters)

expect_true(
  current = all(f1_scores >= 0.95),
  info = "seurat <> bixverse - broadly similar clusters identified"
)
