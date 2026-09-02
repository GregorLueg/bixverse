# edgeR quasi-likelihood tests -------------------------------------------------

source("helper_sc.R", local = TRUE)

## synthetic test data ---------------------------------------------------------

set.seed(42L)

n_features <- 800L
n_samples <- 8L
n_de <- 80L

grp <- rep(c(0, 1), each = n_samples / 2)
base_expr <- 2^runif(n_features, 3, 11)
true_lfc <- c(rnorm(n_de, 1.2, 0.2), rep(0, n_features - n_de))

mu <- outer(base_expr, rep(1, n_samples)) * 2^outer(true_lfc, grp)
counts <- matrix(
  rnbinom(n_features * n_samples, mu = mu, size = 12),
  nrow = n_features
)
rownames(counts) <- sprintf("gene_%04i", seq_len(n_features))
colnames(counts) <- sprintf("sample_%02i", seq_len(n_samples))

design <- stats::model.matrix(~ factor(grp))
colnames(design) <- c("(Intercept)", "trt")

# tests ------------------------------------------------------------------------

## parameters ------------------------------------------------------------------

edger_params <- params_edger_ql()

expect_true(
  current = bixverse:::checkEdgeRQlParams(edger_params),
  info = "the default edgeR params pass their own check"
)

expect_equal(
  current = edger_params$norm_method,
  target = "TMM",
  info = "edgeR params default to TMM"
)

expect_error(
  current = params_edger_ql(norm_method = "nonsense"),
  info = "params_edger_ql rejects an unknown normalisation method"
)

expect_error(
  current = params_edger_ql(min_mean = -1),
  info = "params_edger_ql rejects a negative min_mean"
)

expect_true(
  current = is.character(bixverse:::checkEdgeRQlParams(
    utils::modifyList(edger_params, list(norm_method = "nonsense"))
  )),
  info = "checkEdgeRQlParams returns a message for a bad norm_method"
)

expect_true(
  current = is.character(bixverse:::checkEdgeRQlParams(edger_params[-1])),
  info = "checkEdgeRQlParams returns a message when an element is missing"
)

## input validation ------------------------------------------------------------

expect_error(
  current = run_edger_ql(counts = counts, design = design[, 1, drop = FALSE]),
  pattern = "at least two columns",
  info = "run_edger_ql rejects a design with a single column"
)

expect_error(
  current = run_edger_ql(counts = counts, design = design[1:4, ]),
  pattern = "rows against",
  info = "run_edger_ql rejects a design that does not match the sample count"
)

expect_error(
  current = run_edger_ql(counts = -counts, design = design),
  pattern = "negative values",
  info = "run_edger_ql rejects negative counts"
)

expect_error(
  current = run_edger_ql(
    counts = counts,
    design = design,
    coef = 2L,
    contrast = c(0, 1)
  ),
  pattern = "not both",
  info = "run_edger_ql rejects coef and contrast together"
)

expect_error(
  current = run_edger_ql(counts = counts, design = design, coef = 99L),
  info = "run_edger_ql rejects an out-of-range coefficient"
)

expect_error(
  current = run_edger_ql(
    counts = counts,
    design = design,
    contrast = c(0, 1, 0)
  ),
  pattern = "rows against",
  info = "run_edger_ql rejects a contrast of the wrong length"
)

## results shape ---------------------------------------------------------------

res <- run_edger_ql(counts = counts, design = design)

expect_true(
  current = checkmate::testDataTable(res),
  info = "run_edger_ql returns a data.table"
)

expect_true(
  current = checkmate::testNames(
    names(res),
    permutation.of = c(
      "feature_id",
      "log_fc",
      "log_cpm",
      "f_stat",
      "p_value",
      "fdr"
    )
  ),
  info = "run_edger_ql returns the expected columns"
)

expect_true(
  current = all(res$feature_id %in% rownames(counts)),
  info = "the returned feature ids come from the count matrix rownames"
)

expect_true(
  current = all(res$p_value >= 0 & res$p_value <= 1),
  info = "p-values are in [0, 1]"
)

expect_true(
  current = all(res$fdr >= res$p_value - 1e-12),
  info = "the FDR is never below the raw p-value"
)

## it finds the planted signal -------------------------------------------------

de_genes <- rownames(counts)[seq_len(n_de)]
hits <- res[fdr <= 0.05, feature_id]

# 4 vs 4 samples at a negative binomial size of 12 lands around 0.75 recall
expect_true(
  current = length(intersect(hits, de_genes)) / n_de > 0.7,
  info = "most of the planted DE genes are recovered"
)

expect_true(
  current = length(setdiff(hits, de_genes)) / length(hits) < 0.1,
  info = "few of the calls are false positives"
)

expect_true(
  current = mean(res[feature_id %in% de_genes, log_fc]) > 0.8,
  info = "the log fold changes point the right way"
)

## coefficient selection -------------------------------------------------------

res_by_name <- run_edger_ql(counts = counts, design = design, coef = "trt")
res_by_pos <- run_edger_ql(counts = counts, design = design, coef = 2L)

expect_equal(
  current = res_by_name,
  target = res_by_pos,
  info = "naming the coefficient matches selecting it by position"
)

# `res` picked up an auto-index from the subsetting above, so compare the data
expect_equal(
  current = as.data.frame(res_by_pos),
  target = as.data.frame(res),
  info = "the default coefficient is the last design column"
)

res_contrast <- run_edger_ql(
  counts = counts,
  design = design,
  contrast = c(0, 1)
)

expect_equal(
  current = res_contrast$log_fc,
  target = res$log_fc,
  tolerance = 1e-8,
  info = "a unit contrast on a coefficient reproduces testing that coefficient"
)

## filtering -------------------------------------------------------------------

res_nofilter <- run_edger_ql(
  counts = counts,
  design = design,
  edger_params = params_edger_ql(filter = FALSE)
)

expect_true(
  current = nrow(res_nofilter) >= nrow(res),
  info = "turning filterByExpr off keeps at least as many features"
)

res_minmean <- run_edger_ql(
  counts = counts,
  design = design,
  edger_params = params_edger_ql(filter = FALSE, min_mean = 50)
)

expect_true(
  current = nrow(res_minmean) < nrow(res_nofilter),
  info = "min_mean drops the low-count features"
)

## parity against edgeR --------------------------------------------------------

# The legacy pipeline is the like-for-like comparison: it runs estimateDisp and
# applies the Poisson bound, which is what glmQLFit(legacy = TRUE) does.
edger_ref <- function(counts, design, filter = TRUE) {
  y <- edgeR::DGEList(counts = counts)
  keep <- if (filter) {
    edgeR::filterByExpr(y, design = design)
  } else {
    rep(TRUE, nrow(counts))
  }
  # keep.lib.sizes = TRUE is `[.DGEList`'s default and what the Rust does
  y <- y[keep, , keep.lib.sizes = TRUE]
  y <- edgeR::calcNormFactors(y, method = "TMM")
  y <- edgeR::estimateDisp(y, design)
  fit <- edgeR::glmQLFit(y, design, robust = FALSE, legacy = TRUE)
  list(
    keep = keep,
    res = as.data.frame(edgeR::topTags(
      edgeR::glmQLFTest(fit, coef = ncol(design)),
      sort.by = "none",
      n = Inf
    ))
  )
}

ref <- edger_ref(counts, design, filter = TRUE)
ours <- run_edger_ql(
  counts = counts,
  design = design,
  edger_params = params_edger_ql(filter = TRUE, legacy = TRUE)
)

expect_equal(
  current = ours$feature_id,
  target = rownames(ref$res),
  info = "the Rust filter keeps exactly the features filterByExpr keeps"
)

# Independent implementations of the same chain, so this is agreement rather
# than bit-identity. The log fold changes are closed-form given the fit and
# come out tightest; the p-values inherit the fit's error through the F tail.
expect_equal(
  current = ours$log_fc,
  target = ref$res$logFC,
  tolerance = 1e-6,
  info = "log fold changes match edgeR"
)

expect_equal(
  current = ours$log_cpm,
  target = ref$res$logCPM,
  tolerance = 1e-3,
  info = "average log CPM matches edgeR"
)

expect_equal(
  current = ours$f_stat,
  target = ref$res$F,
  tolerance = 1e-3,
  info = "the quasi-likelihood F statistic matches edgeR"
)

expect_equal(
  current = ours$p_value,
  target = ref$res$PValue,
  tolerance = 1e-3,
  info = "p-values match edgeR"
)

expect_equal(
  current = ours$fdr,
  target = ref$res$FDR,
  tolerance = 1e-3,
  info = "adjusted p-values match edgeR"
)

expect_true(
  current = cor(ours$p_value, ref$res$PValue, method = "spearman") > 0.9999,
  info = "the p-value ranking is the same as edgeR's"
)

## pseudobulk ------------------------------------------------------------------

test_temp_dir <- sc_test_dir("dge_edger")

single_cell_test_data <- sc_test_fixture(
  min_lib_size = 300L,
  min_genes_exp = 45L,
  min_cells_exp = 500L,
  hvg_to_keep = 50L,
  no_pcs = 20L,
  syn_data_params = params_sc_synthetic_data(
    n_samples = 6L,
    sample_bias = "even"
  )
)

sc_object <- sc_test_object(
  test_temp_dir,
  single_cell_test_data,
  sc_qc_param = sc_test_qc_params(single_cell_test_data, target_size = 1000)
)

obs_data <- sc_object[[c("cell_id", "sample_id")]]
cell_list <- split(obs_data$cell_id, obs_data$sample_id)

pb_design <- stats::model.matrix(
  ~ factor(rep(c("ctr", "trt"), each = 3))
)
colnames(pb_design) <- c("(Intercept)", "trt")

pb_res <- pseudobulk_dge_sc(
  object = sc_object,
  cell_list = cell_list,
  design = pb_design,
  edger_params = params_edger_ql(filter = FALSE),
  .verbose = FALSE
)

expect_true(
  current = checkmate::testDataTable(pb_res),
  info = "pseudobulk_dge_sc returns a data.table"
)

expect_true(
  current = all(pb_res$feature_id %in% get_gene_names(sc_object)),
  info = "pseudobulk_dge_sc returns gene identifiers"
)

expect_error(
  current = pseudobulk_dge_sc(
    object = sc_object,
    cell_list = cell_list,
    design = pb_design[1:2, ],
    .verbose = FALSE
  ),
  pattern = "rows against",
  info = "pseudobulk_dge_sc rejects a design that does not match the samples"
)

# the samples carry no real condition effect, so almost nothing should call
expect_true(
  current = sum(pb_res$fdr <= 0.05) < 0.05 * nrow(pb_res),
  info = "few false positives on pseudobulk with no planted effect"
)

# clean up ---------------------------------------------------------------------

sc_test_cleanup(test_temp_dir)
