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
no_pcs <- 20L

## synthetic test data ---------------------------------------------------------

?params_sc_synthetic_data

single_cell_test_data <- generate_single_cell_test_data(
  syn_data_params = params_sc_synthetic_data(
    n_samples = 6L,
    sample_bias = "even"
  )
)

genes_pass <- which(
  Matrix::colSums(single_cell_test_data$counts != 0) >= min_cells_exp
)

cells_pass <- which(
  (Matrix::rowSums(single_cell_test_data$counts[, genes_pass]) >=
    min_lib_size) &
    (Matrix::rowSums(single_cell_test_data$counts[, genes_pass] != 0) >=
      min_genes_exp)
)

## object gen ------------------------------------------------------------------

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

sc_object <- find_neighbours_sc(sc_object)

miloR_obj <- get_miloR_abundances_sc(
  object = sc_object,
  sample_id_col = "sample_id"
)

str(miloR_obj)

design_df <- data.frame(
  grps = c(rep("ctr", 3), rep("trt", 3))
)
rownames(design_df) <- sprintf("sample_%i", 1:6)

x = miloR_obj
design = as.formula(~grps)
design_df
norm_method = c("TMM", "RLE", "logMS")
min_mean = 0
robust = TRUE
fdr_weighting = c("k-distance", "graph-overlap", "none")


norm_method = norm_method[1]
fdr_weighting = fdr_weighting[1]

mm <- stats::model.matrix(design, data = design_df)

if (ncol(x$sample_counts) != nrow(mm)) {
  stop(
    "Design matrix (",
    nrow(mm),
    ") and sample counts (",
    ncol(x$sample_counts),
    ") dimensions don't match"
  )
}

if (any(colnames(x$sample_counts) != rownames(mm))) {
  if (!all(colnames(x$sample_counts) %in% rownames(mm))) {
    stop("Sample names in counts and design matrix don't match")
  }
  warning("Reordering design matrix to match sample counts")
  mm <- mm[colnames(x$sample_counts), ]
}

keep_nh <- if (min_mean > 0) {
  rowMeans(x$sample_counts) >= min_mean
} else {
  rep(TRUE, nrow(x$sample_counts))
}

dge <- edgeR::DGEList(
  counts = x$sample_counts[keep_nh, , drop = FALSE],
  lib.size = colSums(x$sample_counts)
)

if (norm_method %in% c("TMM", "RLE")) {
  dge <- edgeR::calcNormFactors(dge, method = norm_method)
}


dge <- edgeR::estimateDisp(dge, mm)
fit <- edgeR::glmQLFit(dge, mm, robust = robust, legacy = TRUE)


n_coef <- ncol(mm)
res <- edgeR::topTags(
  edgeR::glmQLFTest(fit, coef = n_coef),
  sort.by = "none",
  n = Inf
)
res <- as.data.frame(res)
