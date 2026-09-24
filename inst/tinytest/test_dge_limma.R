# limma / edgeR parity ---------------------------------------------------------

# The bulk DGE workflow runs on edge-rs. edgeR and limma are only Suggests, so
# this file checks the Rust pieces against them when they are installed.

if (
  !requireNamespace("edgeR", quietly = TRUE) ||
    !requireNamespace("limma", quietly = TRUE)
) {
  exit_file("edgeR and limma are not installed")
}

## data ------------------------------------------------------------------------

test_data <- qs2::qs_read("./synthetic_data/dge_test_data.qs")
meta <- test_data$meta_data
counts <- test_data$counts[, meta$sample_id]
group <- meta$dex

# edgeR's defaults
keep_r <- edgeR::filterByExpr(counts, group = group)
counts_kept <- counts[keep_r, ]

## filtering, normalisation, cpm -----------------------------------------------

keep_rs <- rs_filter_by_expr(
  counts = counts,
  group = as.integer(factor(group)),
  lib_size = NULL,
  min_count = 10,
  min_total_count = 15,
  min_prop = 0.7
)

expect_equal(
  current = keep_rs,
  target = unname(keep_r),
  info = "parity - filterByExpr"
)

# edgeR keeps the pre-filter library sizes when subsetting a DGEList
y <- edgeR::DGEList(counts)[keep_r, , keep.lib.sizes = TRUE]
y <- edgeR::normLibSizes(y)

nf_rs <- rs_calc_norm_factors(
  counts = counts_kept,
  lib_size = colSums(counts),
  norm_method = "TMM"
)

expect_equal(
  current = nf_rs,
  target = unname(y$samples$norm.factors),
  tolerance = 1e-10,
  info = "parity - calcNormFactors (TMM)"
)

expect_equal(
  current = rs_cpm(counts, lib_size = NULL, log = TRUE, prior_count = 2),
  target = edgeR::cpm(counts, log = TRUE),
  tolerance = 1e-10,
  info = "parity - cpm (log, keeps dimnames)"
)

## voom ------------------------------------------------------------------------

design <- stats::model.matrix(~ 0 + group)
eff_lib <- y$samples$lib.size * y$samples$norm.factors

voom_rs <- rs_voom_normalise(
  counts = counts_kept,
  design = design,
  lib_size = eff_lib,
  span = 0.5,
  adaptive_span = TRUE
)
voom_r <- limma::voom(y, design, adaptive.span = TRUE)

expect_equal(
  current = voom_rs$e,
  target = voom_r$E,
  tolerance = 1e-10,
  info = "parity - voom E"
)

expect_equal(
  current = unname(voom_rs$weights),
  target = unname(voom_r$weights),
  tolerance = 1e-8,
  info = "parity - voom weights"
)

## limma-voom chain ------------------------------------------------------------

res_rs <- run_limma_voom(
  meta_data = meta,
  main_contrast = "dex",
  counts = counts_kept,
  .verbose = FALSE
)

# same design and contrast run_limma_voom builds: ~ 0 + dex, trt - untrt
y_kept <- edgeR::normLibSizes(edgeR::DGEList(counts_kept))
fit <- edgeR::voomLmFit(y_kept, design, sample.weights = FALSE)
fit <- limma::contrasts.fit(fit, contrasts = c(1, -1))
fit <- limma::eBayes(fit)
res_r <- data.table::as.data.table(
  limma::topTable(fit, number = Inf, sort.by = "none", confint = TRUE),
  keep.rownames = "gene_id"
)

res_rs <- res_rs[match(res_r$gene_id, res_rs$gene_id)]

expect_equal(
  current = res_rs$gene_id,
  target = res_r$gene_id,
  info = "parity - limma-voom genes"
)

for (col in c("logFC", "CI.L", "CI.R", "AveExpr", "t", "B")) {
  expect_equal(
    current = res_rs[[col]],
    target = res_r[[col]],
    tolerance = 1e-8,
    info = sprintf("parity - limma-voom %s", col)
  )
}

for (col in c("P.Value", "adj.P.Val")) {
  expect_equal(
    current = -log10(res_rs[[col]]),
    target = -log10(res_r[[col]]),
    tolerance = 1e-8,
    info = sprintf("parity - limma-voom %s", col)
  )
}

## removeBatchEffect -----------------------------------------------------------

log_cpm <- edgeR::cpm(counts_kept, log = TRUE)
design_bio <- stats::model.matrix(~group)

expect_equal(
  current = rs_remove_batch_effect(
    x = log_cpm,
    batch = as.integer(factor(meta$cell)),
    design = design_bio
  ),
  target = limma::removeBatchEffect(
    log_cpm,
    batch = meta$cell,
    design = design_bio
  ),
  tolerance = 1e-10,
  info = "parity - removeBatchEffect"
)
