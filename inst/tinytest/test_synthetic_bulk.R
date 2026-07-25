# synthetic bulk data generator tests ------------------------------------------

# These guard the generator itself, so the frozen fixtures in the downstream
# module-detection tests do not have to. Thresholds are deliberately loose: they
# should catch a broken generator, not a re-tuned one.

membership <- is_hub <- NULL

## helpers ---------------------------------------------------------------------

# Mean absolute Spearman correlation within modules minus the same across
# modules. This is the number that decides whether any module detection method
# has a chance.
cor_gap <- function(synth) {
  norm_counts <- edgeR::cpm(synth$counts, log = TRUE)
  cor_mat <- stats::cor(t(norm_counts), method = "spearman")
  diag(cor_mat) <- NA

  memb <- synth$module_data$membership
  in_module <- outer(memb, memb, "==") & memb > 0
  across_modules <- outer(memb, memb, "!=") & outer(memb > 0, memb > 0, "&")

  within <- mean(abs(cor_mat[in_module]), na.rm = TRUE)
  across <- mean(abs(cor_mat[across_modules]), na.rm = TRUE)

  list(within = within, across = across, gap = within - across, cor = cor_mat)
}

# Absolute correlation of each module's first principal component against the
# latent factor it was built on.
factor_recovery <- function(synth) {
  norm_counts <- edgeR::cpm(synth$counts, log = TRUE)
  memb <- synth$module_data$membership

  purrr::map_dbl(seq_len(nrow(synth$module_factors)), \(k) {
    module_genes <- which(memb == k)
    pc1 <- stats::prcomp(
      t(norm_counts[module_genes, , drop = FALSE]),
      scale. = TRUE
    )$x[, 1]
    abs(stats::cor(pc1, synth$module_factors[k, ]))
  })
}

excess_kurtosis <- function(x) {
  centred <- x - mean(x)
  mean(centred^4) / (mean(centred^2)^2) - 3
}

## data ------------------------------------------------------------------------

module_sizes <- c(100L, 100L, 100L)

synth_params <- \(generator = "hub_modular", seed = 123L) {
  params_synthetic_bulk_rnaseq(
    num_samples = 100L,
    num_genes = 1000L,
    module_sizes = module_sizes,
    generator = generator,
    seed = seed
  )
}

synth_hub <- synthetic_bulk_cor_matrix(synth_params(generator = "hub_modular"))
synth_modular <- synthetic_bulk_cor_matrix(synth_params(generator = "modular"))
synth_nmf <- synthetic_bulk_cor_matrix(
  synth_params(generator = "non_negative_factor")
)
synth_ica <- synthetic_bulk_cor_matrix(
  synth_params(generator = "non_gaussian_factor")
)

all_generators <- list(
  hub_modular = synth_hub,
  modular = synth_modular,
  non_negative_factor = synth_nmf,
  non_gaussian_factor = synth_ica
)

## object shape ----------------------------------------------------------------

expect_true(
  current = inherits(synth_hub, "synthetic_bulk_data"),
  info = "synthetic bulk - returns a synthetic_bulk_data object"
)

expect_equal(
  current = names(synth_hub),
  target = c("counts", "sparse_counts", "module_data", "module_factors"),
  info = "synthetic bulk - object carries the four expected slots"
)

expect_equal(
  current = names(synth_hub$module_data),
  target = c("gene", "membership", "loading", "is_hub"),
  info = "synthetic bulk - module_data carries the per-gene ground truth"
)

expect_equal(
  current = attr(synth_hub, "synthetic_params")$generator,
  target = "hub_modular",
  info = "synthetic bulk - params are stashed on the object"
)

for (gen in names(all_generators)) {
  synth <- all_generators[[gen]]

  expect_equal(
    current = dim(synth$counts),
    target = c(1000L, 100L),
    info = sprintf("synthetic bulk (%s) - counts are genes x samples", gen)
  )

  expect_true(
    current = all(synth$counts >= 0) &&
      all(synth$counts == round(synth$counts)),
    info = sprintf("synthetic bulk (%s) - counts are non-negative integers", gen)
  )

  expect_equal(
    current = dim(synth$module_factors),
    target = c(length(module_sizes), 100L),
    info = sprintf("synthetic bulk (%s) - factors are modules x samples", gen)
  )

  expect_equal(
    current = colnames(synth$module_factors),
    target = colnames(synth$counts),
    info = sprintf("synthetic bulk (%s) - factor columns match samples", gen)
  )
}

## module ground truth ---------------------------------------------------------

expect_equal(
  current = synth_hub$module_data[
    membership > 0,
    .N,
    by = membership
  ][order(membership)]$N,
  target = module_sizes,
  info = "synthetic bulk - module block sizes match module_sizes"
)

expect_equal(
  current = synth_hub$module_data[membership == 0, .N],
  target = 1000L - sum(module_sizes),
  info = "synthetic bulk - the remainder is background"
)

expect_true(
  current = all(
    synth_hub$module_data[membership > 0, gene] ==
      sprintf("gene_%i", seq_len(sum(module_sizes)))
  ),
  info = "synthetic bulk - modules occupy contiguous blocks from the first gene"
)

expect_true(
  current = all(synth_hub$module_data[membership == 0, loading] == 0) &&
    !any(synth_hub$module_data[membership == 0, is_hub]),
  info = "synthetic bulk - background genes have no loading and are never hubs"
)

expect_true(
  current = all(synth_hub$module_data[membership > 0, loading] > 0),
  info = "synthetic bulk - module genes all carry a positive loading"
)

## correlation structure -------------------------------------------------------

for (gen in names(all_generators)) {
  gap_res <- cor_gap(all_generators[[gen]])

  expect_true(
    current = gap_res$gap > 0.15,
    info = sprintf(
      "synthetic bulk (%s) - within minus cross module correlation gap > 0.15 (got %.3f)",
      gen,
      gap_res$gap
    )
  )
}

## factor recovery -------------------------------------------------------------

for (gen in names(all_generators)) {
  recovery <- factor_recovery(all_generators[[gen]])

  expect_true(
    current = all(recovery > 0.7),
    info = sprintf(
      "synthetic bulk (%s) - module PC1 recovers the latent factor (min %.2f)",
      gen,
      min(recovery)
    )
  )
}

## hubs ------------------------------------------------------------------------

# Hubs are ranked on loading globally across all module genes, not per module,
# so the per-module counts do not have to be equal.
expect_equal(
  current = synth_hub$module_data[, sum(is_hub)],
  target = ceiling(sum(module_sizes) * 0.1),
  info = "synthetic bulk - hub count follows hub_percentile"
)

hub_gap <- cor_gap(synth_hub)
hub_adjacency <- abs(hub_gap$cor) >= 0.5
hub_degree <- rowSums(hub_adjacency, na.rm = TRUE)
hub_flag <- synth_hub$module_data$is_hub
in_module <- synth_hub$module_data$membership > 0

expect_true(
  current = mean(hub_degree[hub_flag]) >
    mean(hub_degree[in_module & !hub_flag]),
  info = "synthetic bulk (hub_modular) - hubs are more connected than non-hubs"
)

expect_equal(
  current = synth_modular$module_data[, sum(is_hub)],
  target = 0L,
  info = "synthetic bulk (modular) - plants no hubs"
)

## generator-specific factor properties ----------------------------------------

expect_true(
  current = all(synth_nmf$module_factors >= 0),
  info = "synthetic bulk (non_negative_factor) - factors are non-negative"
)

expect_true(
  current = excess_kurtosis(as.vector(synth_ica$module_factors)) > 1,
  info = "synthetic bulk (non_gaussian_factor) - factors are heavy tailed"
)

expect_true(
  current = abs(excess_kurtosis(as.vector(synth_hub$module_factors))) < 1,
  info = "synthetic bulk (hub_modular) - factors are approximately Gaussian"
)

## reproducibility -------------------------------------------------------------

expect_identical(
  current = synthetic_bulk_cor_matrix(synth_params(seed = 7L))$counts,
  target = synthetic_bulk_cor_matrix(synth_params(seed = 7L))$counts,
  info = "synthetic bulk - the same seed reproduces the same counts"
)

expect_false(
  current = identical(
    synthetic_bulk_cor_matrix(synth_params(seed = 7L))$counts,
    synthetic_bulk_cor_matrix(synth_params(seed = 8L))$counts
  ),
  info = "synthetic bulk - a different seed gives different counts"
)

## no modules ------------------------------------------------------------------

synth_none <- synthetic_bulk_cor_matrix(
  params_synthetic_bulk_rnaseq(num_genes = 200L, module_sizes = integer(0))
)

expect_true(
  current = all(synth_none$module_data$membership == 0) &&
    all(synth_none$module_data$loading == 0) &&
    !any(synth_none$module_data$is_hub),
  info = "synthetic bulk - an empty module_sizes plants no structure"
)

expect_equal(
  current = dim(synth_none$module_factors),
  target = c(0L, 100L),
  info = "synthetic bulk - no modules gives a zero-row factor matrix"
)

## parameter validation --------------------------------------------------------

# `module_sizes` reaches Rust through `as_integer_slice()`, where a double
# vector reads as absent and silently falls back to the crate default. Likewise
# an unknown `generator` silently resolves to "hub_modular". Both must be caught
# in R.
expect_error(
  current = params_synthetic_bulk_rnaseq(module_sizes = c(100, 100, 100)),
  info = "synthetic bulk params - a double module_sizes is rejected"
)

expect_error(
  current = params_synthetic_bulk_rnaseq(module_sizes = NULL),
  info = "synthetic bulk params - a NULL module_sizes is rejected"
)

expect_error(
  current = params_synthetic_bulk_rnaseq(generator = "typo"),
  info = "synthetic bulk params - an unknown generator is rejected"
)

expect_error(
  current = params_synthetic_bulk_rnaseq(
    num_genes = 100L,
    module_sizes = c(60L, 60L)
  ),
  info = "synthetic bulk params - module sizes exceeding num_genes are rejected"
)

expect_error(
  current = params_synthetic_bulk_rnaseq(seed = -1L),
  info = "synthetic bulk params - a negative seed is rejected"
)

expect_error(
  current = params_synthetic_bulk_rnaseq(hub_percentile = 0),
  info = "synthetic bulk params - a zero hub_percentile is rejected"
)

expect_error(
  current = assertSyntheticBulkParams(list(nonsense = 1)),
  info = "synthetic bulk params - assertion rejects a list missing the names"
)

expect_error(
  current = assertSyntheticBulkParams(
    utils::modifyList(
      params_synthetic_bulk_rnaseq(),
      list(hub_percentile = 2.0)
    )
  ),
  info = "synthetic bulk params - assertion rejects an out-of-range value"
)

expect_error(
  current = params_bulk_sparsity(target_library_size = -1),
  info = "bulk sparsity params - a negative library size is rejected"
)

expect_error(
  current = assertBulkSparsityParams(list(strategy = "seq_depth")),
  info = "bulk sparsity params - assertion rejects a list missing the names"
)

## dropout ---------------------------------------------------------------------

sparse_synth <- simulate_dropouts(synth_hub)

expect_equal(
  current = dim(sparse_synth$sparse_counts),
  target = dim(synth_hub$counts),
  info = "dropout - the sparsified matrix keeps its shape"
)

expect_identical(
  current = dimnames(sparse_synth$sparse_counts),
  target = dimnames(synth_hub$counts),
  info = "dropout - the sparsified matrix keeps its dimnames"
)

expect_true(
  current = all(
    colSums(sparse_synth$sparse_counts) <= colSums(sparse_synth$counts)
  ),
  info = "dropout - library sizes never grow"
)

expect_true(
  current = mean(sparse_synth$sparse_counts == 0) >
    mean(sparse_synth$counts == 0),
  info = "dropout - sparsity increases"
)

expect_identical(
  current = simulate_dropouts(synth_hub)$sparse_counts,
  target = sparse_synth$sparse_counts,
  info = "dropout - the same seed reproduces the same matrix"
)

expect_false(
  current = identical(
    simulate_dropouts(synth_hub, params_bulk_sparsity(seed = 99L))$sparse_counts,
    sparse_synth$sparse_counts
  ),
  info = "dropout - a different seed gives a different matrix"
)

expect_true(
  current = calculate_sparsity_stats(sparse_synth)$added_sparsity > 0,
  info = "dropout - calculate_sparsity_stats reports added sparsity"
)

# Dropout should thin the counts without destroying the planted structure.
sparse_gap <- cor_gap(list(
  counts = sparse_synth$sparse_counts,
  module_data = sparse_synth$module_data
))

expect_true(
  current = sparse_gap$gap > 0.1,
  info = sprintf(
    "dropout - the module correlation structure survives (gap %.3f)",
    sparse_gap$gap
  )
)
