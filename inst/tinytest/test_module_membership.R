# module membership from factorisation loadings --------------------------------

# Guards the shared helper behind get_modules() for ICA, NMF and DGRDL. The
# property that matters is that membership is sparse and NOT exclusive: a gene
# loading on several components belongs to several modules, and a gene in no tail
# belongs to none. An argmax assignment gets both of those wrong.

gene <- module_id <- N <- NULL

## params ----------------------------------------------------------------------

expect_equal(
  current = params_module_membership(),
  target = list(method = "zscore", cutoff = 3.0, fdr = 0.05, tails = "auto"),
  info = "module membership params - defaults"
)

expect_error(
  current = params_module_membership(method = "typo"),
  info = "module membership params - unknown method rejected"
)

expect_error(
  current = params_module_membership(tails = "sideways"),
  info = "module membership params - unknown tails rejected"
)

expect_error(
  current = params_module_membership(cutoff = 0),
  info = "module membership params - non-positive cutoff rejected"
)

expect_error(
  current = params_module_membership(fdr = 1.5),
  info = "module membership params - fdr above 1 rejected"
)

expect_error(
  current = assertModuleMembershipParams(list(method = "zscore")),
  info = "module membership params - assertion rejects an incomplete list"
)

## synthetic loading matrix ----------------------------------------------------

# 200 genes, 3 components. Genes 1:20 load hard on component 1, 15:34 on
# component 2 (so 15:20 sit in BOTH), 40:59 load hard and negative on
# component 3. Everything else is noise.
set.seed(42L)
loadings <- matrix(
  rnorm(200 * 3, sd = 0.1),
  nrow = 200,
  dimnames = list(sprintf("gene_%i", 1:200), sprintf("comp_%i", 1:3))
)
loadings[1:20, 1] <- 5
loadings[15:34, 2] <- 5
loadings[40:59, 3] <- -5

mods <- bixverse:::.modules_from_loadings(loadings)

expect_true(
  current = data.table::is.data.table(mods),
  info = "membership - returns a data.table"
)

expect_equal(
  current = names(mods),
  target = c("gene", "module_id", "loading", "sign", "z"),
  info = "membership - zscore method reports the z score"
)

expect_true(
  current = mods[, .N, by = gene][N > 1, .N] > 0,
  info = "membership - a gene can belong to more than one module"
)

expect_true(
  current = all(sprintf("gene_%i", 15:20) %in%
    mods[module_id == "comp_1", gene]) &&
    all(sprintf("gene_%i", 15:20) %in% mods[module_id == "comp_2", gene]),
  info = "membership - the overlapping block lands in both components"
)

expect_true(
  current = data.table::uniqueN(mods$gene) < nrow(loadings),
  info = "membership - noise genes are excluded, giving a real background"
)

expect_true(
  current = all(mods[module_id == "comp_3", sign] == "neg"),
  info = "membership - the negative block is signed neg"
)

# Only the planted block is asserted: under a two-sided rule, noise genes sitting
# in comp_1's lower tail legitimately pass and are signed neg.
expect_true(
  current = all(
    mods[module_id == "comp_1" & gene %in% sprintf("gene_%i", 1:20), sign] ==
      "pos"
  ),
  info = "membership - the planted positive block is signed pos"
)

# Ordered by component then descending absolute loading.
expect_true(
  current = mods[module_id == "comp_1", all(diff(abs(loading)) <= 0)],
  info = "membership - rows are ordered by descending absolute loading"
)

## cutoff behaviour ------------------------------------------------------------

strict <- bixverse:::.modules_from_loadings(
  loadings,
  params_module_membership(cutoff = 10)
)
loose <- bixverse:::.modules_from_loadings(
  loadings,
  params_module_membership(cutoff = 1)
)

expect_true(
  current = nrow(strict) < nrow(mods) && nrow(mods) < nrow(loose),
  info = "membership - a higher cutoff keeps fewer genes"
)

## fdr method ------------------------------------------------------------------

mods_fdr <- bixverse:::.modules_from_loadings(
  loadings,
  params_module_membership(method = "fdr")
)

expect_equal(
  current = names(mods_fdr),
  target = c("gene", "module_id", "loading", "sign", "padj"),
  info = "membership - fdr method reports the adjusted p-value"
)

expect_true(
  current = all(mods_fdr$padj < 0.05),
  info = "membership - fdr method respects its threshold"
)

expect_true(
  current = all(sprintf("gene_%i", 1:14) %in%
    mods_fdr[module_id == "comp_1", gene]),
  info = "membership - fdr method recovers the planted block"
)

## tails -----------------------------------------------------------------------

non_negative <- abs(loadings)

expect_true(
  current = all(
    bixverse:::.modules_from_loadings(non_negative)$sign == "pos"
  ),
  info = "membership - auto tails goes one-sided on non-negative loadings"
)

expect_true(
  current = nrow(bixverse:::.modules_from_loadings(
    loadings,
    params_module_membership(tails = "upper")
  )) <
    nrow(bixverse:::.modules_from_loadings(
      loadings,
      params_module_membership(tails = "both")
    )),
  info = "membership - forcing the upper tail drops the negative block"
)

## degenerate input ------------------------------------------------------------

flat <- matrix(
  1,
  nrow = 10,
  ncol = 2,
  dimnames = list(sprintf("gene_%i", 1:10), c("comp_1", "comp_2"))
)

expect_warning(
  current = bixverse:::.modules_from_loadings(flat),
  info = "membership - a component with no spread warns rather than erroring"
)

expect_equal(
  current = nrow(suppressWarnings(bixverse:::.modules_from_loadings(flat))),
  target = 0L,
  info = "membership - a degenerate matrix yields no modules"
)

expect_error(
  current = bixverse:::.modules_from_loadings(unname(loadings)),
  info = "membership - a loading matrix without rownames is rejected"
)
