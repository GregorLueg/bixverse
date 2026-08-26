# bulk NMF tests ---------------------------------------------------------------

library(magrittr)

## data ------------------------------------------------------------------------

# The Gamma-factor generator gives NMF a non-negative ground truth it can
# actually reach. All four generators land within 0.71 to 0.83 recall here, so
# this is the theoretically matched choice rather than the flattering one.
synth <- synthetic_bulk_cor_matrix(
  params_synthetic_bulk_rnaseq(
    num_samples = 60L,
    num_genes = 300L,
    module_sizes = c(50L, 50L, 50L),
    generator = "non_negative_factor",
    seed = 123L
  )
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
  current = all(
    c("gene", "module_id", "loading", "sign") %in%
      names(nmf_modules_dt)
  ),
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
  current = all(
    c("losses", "converged", "best_idx") %in%
      names(nmf_stab_diagnostics)
  ),
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
  current = recall >= 0.6,
  info = sprintf(
    "nmf_bulk - recall vs planted modules >= 0.6 (got %.3f)",
    recall
  )
)

## membership is sparse and non-exclusive --------------------------------------

expect_true(
  current = data.table::uniqueN(nmf_modules_dt$gene) < nrow(truth),
  info = "nmf_bulk - membership is sparse, background genes are not assigned"
)

expect_true(
  current = all(nmf_modules_dt$sign == "pos"),
  info = "nmf_bulk - non-negative loadings give an upper-tail-only test"
)

expect_true(
  current = "z" %in% names(nmf_modules_dt),
  info = "nmf_bulk - the zscore threshold statistic is reported"
)

# A tighter cutoff must keep strictly fewer genes.
nmf_strict <- nmf_bulk(
  BulkCoExp(raw_data = t(counts), meta_data = meta_data) %>%
    preprocess_bulk_coexp(scaling = FALSE, hvg = NULL, .verbose = FALSE),
  k = 5L,
  membership_params = params_module_membership(cutoff = 6),
  seed = 42L,
  .verbose = FALSE
)

expect_true(
  current = nrow(get_modules(get_results(nmf_strict))) < nrow(nmf_modules_dt),
  info = "nmf_bulk - a stricter membership cutoff keeps fewer genes"
)

## consensus nmf ---------------------------------------------------------------

# The density filter is off for the tests. resolve_n_neighbours is
# ceiling(0.3 * n_runs), so at a small n_runs a live filter is jumpy enough on
# a noisy fixture to drop below k survivors, which is a hard error and a flaky
# test. The filter itself gets its own deliberate test below.
consensus_params <- params_nmf_consensus(density_threshold = 2)

nmf_consensus_obj <- BulkCoExp(
  raw_data = t(counts),
  meta_data = meta_data
) %>%
  preprocess_bulk_coexp(scaling = FALSE, hvg = NULL, .verbose = FALSE) %>%
  consensus_nmf_bulk(
    k = 5L,
    n_runs = 6L,
    nmf_consensus_params = consensus_params,
    seed = 42L,
    .verbose = FALSE
  )

nmf_consensus_result <- get_results(nmf_consensus_obj)

expect_true(
  current = inherits(nmf_consensus_result, "BulkModuleResult"),
  info = "consensus_nmf_bulk - final_results is a BulkModuleResult"
)

consensus_gene_loadings <- get_factors(
  nmf_consensus_result,
  which = "gene_loadings"
)
consensus_sample_activity <- get_factors(
  nmf_consensus_result,
  which = "sample_activity"
)

expect_equal(
  current = dim(consensus_gene_loadings),
  target = c(nrow(counts), 5L),
  info = "consensus_nmf_bulk - gene_loadings has (n_genes, k) shape"
)

expect_equal(
  current = dim(consensus_sample_activity),
  target = c(ncol(counts), 5L),
  info = "consensus_nmf_bulk - sample_activity has (n_samples, k) shape"
)

expect_true(
  current = all(consensus_gene_loadings >= 0) &&
    all(consensus_sample_activity >= 0),
  info = "consensus_nmf_bulk - both factors are non-negative"
)

consensus_stability <- get_nmf_stability(nmf_consensus_obj)

expect_true(
  current = all(
    c(
      "stability",
      "rel_error",
      "rel_run_errors",
      "clusters",
      "cluster_sizes",
      "n_dropped",
      "n_empty_clusters"
    ) %in%
      names(consensus_stability)
  ),
  info = "consensus_nmf_bulk - get_nmf_stability returns consensus diagnostics"
)

expect_true(
  current = consensus_stability$stability >= -1 &&
    consensus_stability$stability <= 1,
  info = "consensus_nmf_bulk - stability is a silhouette in [-1, 1]"
)

expect_equal(
  current = length(consensus_stability$rel_run_errors),
  target = 6L,
  info = "consensus_nmf_bulk - one relative error per restart"
)

expect_equal(
  current = nrow(consensus_stability$clusters),
  target = 5L * 6L,
  info = "consensus_nmf_bulk - one cluster row per pooled component"
)

expect_true(
  current = all(
    c(
      "component_id",
      "run",
      "component",
      "pooled_idx",
      "cluster",
      "local_density",
      "silhouette",
      "kept"
    ) %in%
      names(consensus_stability$clusters)
  ),
  info = "consensus_nmf_bulk - cluster table carries the per-component fields"
)

expect_equal(
  current = sum(consensus_stability$cluster_sizes$n),
  target = sum(consensus_stability$clusters$kept),
  info = "consensus_nmf_bulk - cluster sizes account for every survivor"
)

# The stabilised getter must keep working on a stabilised fit, i.e. the
# consensus branch has not swallowed it.
expect_equal(
  current = names(get_nmf_stability(nmf_stab_test)),
  target = c("losses", "converged", "best_idx"),
  info = "get_nmf_stability - stabilised fits still return the old fields"
)

# Same seed, same answer.
nmf_consensus_repeat <- BulkCoExp(
  raw_data = t(counts),
  meta_data = meta_data
) %>%
  preprocess_bulk_coexp(scaling = FALSE, hvg = NULL, .verbose = FALSE) %>%
  consensus_nmf_bulk(
    k = 5L,
    n_runs = 6L,
    nmf_consensus_params = consensus_params,
    seed = 42L,
    .verbose = FALSE
  )

expect_equal(
  current = get_factors(get_results(nmf_consensus_repeat), "gene_loadings"),
  target = consensus_gene_loadings,
  info = "consensus_nmf_bulk - the same seed reproduces the same loadings"
)

consensus_recall <- mean(best_jaccard(get_modules(nmf_consensus_result)))

expect_true(
  current = consensus_recall >= 0.6,
  info = sprintf(
    "consensus_nmf_bulk - recall vs planted modules >= 0.6 (got %.3f)",
    consensus_recall
  )
)

## consensus nmf validation ----------------------------------------------------

expect_error(
  current = consensus_nmf_bulk(nmf_test, k = 1L, n_runs = 4L, .verbose = FALSE),
  info = "consensus_nmf_bulk - k below 2 is refused"
)

expect_error(
  current = consensus_nmf_bulk(nmf_test, k = 3L, n_runs = 1L, .verbose = FALSE),
  info = "consensus_nmf_bulk - fewer than 2 restarts is refused"
)

expect_error(
  current = params_nmf_consensus(density_threshold = 3),
  info = "params_nmf_consensus - a density threshold above 2 is refused"
)

expect_error(
  current = params_nmf_consensus(consensus_target = "garbage"),
  info = "params_nmf_consensus - an unknown consensus target is refused"
)

expect_error(
  current = assertNmfConsensus(list(consensus_target = "h")),
  info = "assertNmfConsensus - an incomplete parameter list is refused"
)

# A filter this tight cannot leave k components standing, and that has to be a
# hard error rather than a quietly degenerate fit.
expect_error(
  current = consensus_nmf_bulk(
    nmf_test,
    k = 5L,
    n_runs = 4L,
    nmf_consensus_params = params_nmf_consensus(density_threshold = 1e-9),
    .verbose = FALSE
  ),
  pattern = "density",
  info = "consensus_nmf_bulk - too tight a density filter errors informatively"
)

## k sweep ---------------------------------------------------------------------

k_sweep <- nmf_k_sweep_bulk(
  nmf_test,
  k_range = 2:5,
  n_runs = 4L,
  nmf_consensus_params = consensus_params,
  seed = 42L,
  .verbose = FALSE
)

expect_true(
  current = inherits(k_sweep, "NmfKSweepResult") &&
    data.table::is.data.table(k_sweep),
  info = "nmf_k_sweep_bulk - returns a data.table-backed NmfKSweepResult"
)

expect_equal(
  current = k_sweep$k,
  target = 2:5,
  info = "nmf_k_sweep_bulk - one row per k, in the order requested"
)

expect_true(
  current = all(
    c(
      "stability",
      "best_error",
      "median_error",
      "consensus_failed",
      "n_dropped",
      "n_empty_clusters",
      "n_converged"
    ) %in%
      names(k_sweep)
  ),
  info = "nmf_k_sweep_bulk - the diagnostic columns are all present"
)

expect_true(
  current = all(k_sweep$best_error <= k_sweep$median_error),
  info = "nmf_k_sweep_bulk - the best restart is never worse than the median"
)

expect_true(
  current = all(diff(k_sweep$best_error) <= 0),
  info = "nmf_k_sweep_bulk - reconstruction error falls as k grows"
)

expect_true(
  current = inherits(plot(k_sweep), "ggplot"),
  info = "nmf_k_sweep_bulk - plot() gives a ggplot"
)

expect_error(
  current = nmf_k_sweep_bulk(nmf_test, k_range = 1:3, .verbose = FALSE),
  info = "nmf_k_sweep_bulk - a k range reaching below 2 is refused"
)

expect_error(
  current = nmf_k_sweep_bulk(nmf_test, k_range = integer(), .verbose = FALSE),
  info = "nmf_k_sweep_bulk - an empty k range is refused"
)
