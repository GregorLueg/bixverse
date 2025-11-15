# differential abundance tests -------------------------------------------------

library(scuttle)
library(scran)
library(scater)


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

# expect_true(
#   current = length(genes_pass) > 80 & length(genes_pass) != 100,
#   info = "sc processing - sensible amount of genes pass"
# )

# expect_true(
#   current = length(cells_pass) > 800 & length(cells_pass) != 1000,
#   info = "sc processing - sensible amount of cells pass"
# )

### bixverse -------------------------------------------------------------------

sc_object <- single_cell_exp(dir_data = tempdir())

sc_object <- load_r_data(
  object = sc_object,
  counts = single_cell_test_data$counts,
  obs = single_cell_test_data$obs,
  var = single_cell_test_data$var,
  sc_qc_param = params_sc_min_quality(
    min_unique_genes = min_genes_exp,
    min_lib_size = min_lib_size,
    min_cells = min_cells_exp,
    target_size = 1000
  ),
  streaming = FALSE,
  .verbose = FALSE
)

sc_object <- find_hvg_sc(
  object = sc_object,
  hvg_no = hvg_to_keep,
  .verbose = FALSE
)


sc_object <- calculate_pca_sc(
  sc_object,
  no_pcs = no_pcs,
  .verbose = FALSE,
  randomised_svd = TRUE
)

### SCE ------------------------------------------------------------------------

sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = Matrix::t(single_cell_test_data$counts)),
  colData = single_cell_test_data$obs,
  rowData = single_cell_test_data$var
)

n_cells_expressing <- rowSums(counts(sce) > 0)

genes_to_keep <- n_cells_expressing >= min_cells_exp

sce <- sce[genes_to_keep, ]

lib_sizes <- colSums(counts(sce))
n_features <- colSums(counts(sce) > 0)

cells_to_keep <- lib_sizes >= min_lib_size & n_features >= min_genes_exp

sce <- sce[, cells_to_keep]

lib_sizes <- colSums(counts(sce))

size_factors <- lib_sizes / 10000
sizeFactors(sce) <- size_factors
sce <- logNormCounts(sce)

gene_var <- modelGeneVar(sce)
hvgs <- getTopHVGs(gene_var, n = hvg_to_keep)

sce <- runPCA(
  sce,
  subset_row = hvgs,
  ncomponents = no_pcs,
  scale = TRUE
)

plot(
  sce@assays@data$logcounts[, 3],
  as.numeric(sc_object[3L, , assay = "norm"])
)


cor(get_pca_factors(sc_object), reducedDim(sce, "PCA"))
