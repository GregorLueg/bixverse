# bulk NMF tests ---------------------------------------------------------------

library(magrittr)

## data ------------------------------------------------------------------------

synth <- synthetic_bulk_cor_matrix(
  num_samples = 60L,
  num_genes = 300L,
  add_modules = TRUE,
  module_sizes = c(50L, 50L, 50L),
  seed = 123L
)

counts <- synth$counts
truth <- synth$module_data

meta_data <- data.table::data.table(sample_id = colnames(counts))

## single-run nmf_bulk ---------------------------------------------------------

nmf_test <- BulkCoExp(
  raw_data = t(counts),
  meta_data = meta_data
) %>%
  preprocess_bulk_coexp(
    scaling = FALSE,
    hvg = NULL,
    .verbose = FALSE
  ) %>%
  nmf_bulk(k = 5L, seed = 42L, .verbose = FALSE)

nmf_result <- get_results(nmf_test)

expect_true(
  current = inherits(nmf_result, "BulkModuleResult"),
  info = "nmf_bulk - final_results is a BulkModuleResult"
)

nmf_modules_dt <- get_modules(nmf_result)

expect_true(
  current = data.table::is.data.table(nmf_modules_dt),
  info = "nmf_bulk - get_modules() returns a data.table"
)

expect_true(
  current = all(c("gene", "module_id", "loading", "sign") %in%
    names(nmf_modules_dt)),
  info = "nmf_bulk - modules DT has gene/module_id/loading/sign"
)

nmf_gene_loadings <- get_factors(nmf_result, which = "gene_loadings")
nmf_sample_activity <- get_factors(nmf_result, which = "sample_activity")

expect_true(
  current = all(nmf_gene_loadings >= 0),
  info = "nmf_bulk - gene_loadings are non-negative"
)

expect_true(
  current = all(nmf_sample_activity >= 0),
  info = "nmf_bulk - sample_activity is non-negative"
)

expect_equal(
  current = dim(nmf_gene_loadings),
  target = c(nrow(counts), 5L),
  info = "nmf_bulk - gene_loadings has (n_genes, k) shape"
)

expect_equal(
  current = dim(nmf_sample_activity),
  target = c(ncol(counts), 5L),
  info = "nmf_bulk - sample_activity has (n_samples, k) shape"
)

## legacy shims ----------------------------------------------------------------

expect_equal(
  current = get_nmf_gene_loadings(nmf_test),
  target = nmf_gene_loadings,
  info = "nmf_bulk - get_nmf_gene_loadings shim equivalent to get_factors"
)

expect_equal(
  current = get_nmf_sample_activity(nmf_test),
  target = nmf_sample_activity,
  info = "nmf_bulk - get_nmf_sample_activity shim equivalent to get_factors"
)

expect_equal(
  current = get_nmf_modules(nmf_test),
  target = nmf_modules_dt,
  info = "nmf_bulk - get_nmf_modules shim equivalent to get_modules"
)

## stabilised_nmf_bulk ---------------------------------------------------------

nmf_stab_test <- BulkCoExp(
  raw_data = t(counts),
  meta_data = meta_data
) %>%
  preprocess_bulk_coexp(
    scaling = FALSE,
    hvg = NULL,
    .verbose = FALSE
  ) %>%
  stabilised_nmf_bulk(
    k = 5L,
    n_runs = 3L,
    seed = 42L,
    .verbose = FALSE
  )

nmf_stab_result <- get_results(nmf_stab_test)

expect_true(
  current = inherits(nmf_stab_result, "BulkModuleResult"),
  info = "stabilised_nmf_bulk - final_results is a BulkModuleResult"
)

nmf_stab_diagnostics <- get_diagnostics(nmf_stab_result)

expect_true(
  current = all(c("losses", "converged", "best_idx") %in%
    names(nmf_stab_diagnostics)),
  info = "stabilised_nmf_bulk - diagnostics has losses/converged/best_idx"
)

expect_equal(
  current = length(nmf_stab_diagnostics$losses),
  target = 3L,
  info = "stabilised_nmf_bulk - one loss per run"
)

## recall vs planted modules ---------------------------------------------------

truth_split <- split(truth$gene, truth$membership)
truth_split <- truth_split[names(truth_split) != "0"]

best_jaccard <- function(assignments) {
  detected <- split(assignments$gene, assignments$module_id)
  purrr::map_dbl(truth_split, \(planted) {
    max(purrr::map_dbl(detected, \(det) {
      length(intersect(planted, det)) / length(union(planted, det))
    }))
  })
}

recall <- mean(best_jaccard(nmf_modules_dt))

expect_true(
  current = recall >= 0.5,
  info = sprintf(
    "nmf_bulk - recall vs planted modules >= 0.5 (got %.3f)",
    recall
  )
)
