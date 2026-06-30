# singscore tests --------------------------------------------------------------

## data ------------------------------------------------------------------------

set.seed(42L)

no_genes <- 200L
no_samples <- 10L
no_gs <- 5L

data <- matrix(
  data = rnorm(no_genes * no_samples, sd = 2),
  nrow = no_genes,
  ncol = no_samples
)

rownames(data) <- sprintf("gene_%i", 1:no_genes)
colnames(data) <- sprintf("sample_%i", 1:no_samples)

# Single signature
up_set <- sample(rownames(data), 20)
down_set <- sample(setdiff(rownames(data), up_set), 20)

# Multi signatures
up_pathways <- purrr::map(1:no_gs, ~ sample(rownames(data), 15))
names(up_pathways) <- sprintf("gs_%i", 1:no_gs)

down_pathways <- purrr::map(1:no_gs, ~ sample(rownames(data), 15))
names(down_pathways) <- sprintf("gs_%i", 1:no_gs)

# Stable genes
stable_genes <- sample(rownames(data), 25)

# Non-matching set
bad_set <- sprintf("invalid_gene_%i", 1:10)

## tests ranking ---------------------------------------------------------------

### rust -----------------------------------------------------------------------

ranks_rs <- rs_rank_matrix_col(exp = data)

expect_equal(
  current = dim(ranks_rs),
  target = dim(data),
  info = "singscore - rank shape preserved"
)

expect_true(
  current = all(ranks_rs >= 1 & ranks_rs <= no_genes),
  info = "singscore - ranks in [1, n_genes]"
)

stable_idx <- match(stable_genes, rownames(data))
stable_ranks_rs <- rs_rank_matrix_col_stable(
  exp = data,
  stable_gene_indices = as.integer(stable_idx)
)

expect_equal(
  current = dim(stable_ranks_rs),
  target = dim(data),
  info = "singscore - stable rank shape preserved"
)

expect_true(
  current = all(stable_ranks_rs > 0 & stable_ranks_rs <= 1),
  info = "singscore - stable ranks normalised to (0, 1]"
)

### r wrapper ------------------------------------------------------------------

ranks_r <- calc_singscore_rank(exp = data)

expect_false(
  current = attr(ranks_r, "stable"),
  info = "singscore - standard rank attribute"
)

expect_equal(
  current = rownames(ranks_r),
  target = rownames(data),
  info = "singscore - rownames preserved"
)

stable_ranks_r <- calc_singscore_rank(exp = data, stable_genes = stable_genes)

expect_true(
  current = attr(stable_ranks_r, "stable"),
  info = "singscore - stable rank attribute"
)

## tests single gene set -------------------------------------------------------

ranks <- calc_singscore_rank(exp = data)

### up only --------------------------------------------------------------------

score_up <- calc_singscore(ranks = ranks, up_set = up_set)

expect_true(
  current = inherits(score_up, "data.table"),
  info = "singscore single - returns data.table"
)

expect_equal(
  current = nrow(score_up),
  target = no_samples,
  info = "singscore single - one row per sample"
)

expect_true(
  current = all(
    c("total_score", "total_dispersion", "sample_id") %in% names(score_up)
  ),
  info = "singscore single - core columns present"
)

expect_false(
  current = any(c("up_score", "down_score") %in% names(score_up)),
  info = "singscore single - no up/down columns when down_set is NULL"
)

expect_true(
  current = all(score_up$total_score >= -0.5 & score_up$total_score <= 0.5),
  info = "singscore single - centered up-only scores in [-0.5, 0.5]"
)

### up + down ------------------------------------------------------------------

score_updn <- calc_singscore(
  ranks = ranks,
  up_set = up_set,
  down_set = down_set
)

expect_true(
  current = all(
    c(
      "total_score",
      "total_dispersion",
      "up_score",
      "up_dispersion",
      "down_score",
      "down_dispersion"
    ) %in%
      names(score_updn)
  ),
  info = "singscore single - all 6 columns with down_set"
)

# Internal consistency: total = up + down
expect_equivalent(
  current = score_updn$total_score,
  target = score_updn$up_score + score_updn$down_score,
  info = "singscore single - total_score = up_score + down_score",
  tolerance = 1e-12
)

# Internal consistency: total dispersion = mean of up and down
expect_equivalent(
  current = score_updn$total_dispersion,
  target = (score_updn$up_dispersion + score_updn$down_dispersion) / 2,
  info = "singscore single - total_dispersion = mean(up, down)",
  tolerance = 1e-12
)

### permutation test ----------------------------------------------------------

n_perm <- 100L

score_perm <- calc_singscore(
  ranks = ranks,
  up_set = up_set,
  down_set = down_set,
  n_permutations = n_perm,
  seed = 42L
)

expect_true(
  current = "pval" %in% names(score_perm),
  info = "singscore single - pval column added with permutations"
)

expect_true(
  current = all(score_perm$pval >= 1 / n_perm & score_perm$pval <= 1),
  info = "singscore single - p-values in [1/B, 1]"
)

expect_equal(
  current = dim(attr(score_perm, "null_distribution")),
  target = c(n_perm, no_samples),
  info = "singscore single - null distribution shape"
)

# Reproducibility
score_perm_2 <- calc_singscore(
  ranks = ranks,
  up_set = up_set,
  down_set = down_set,
  n_permutations = n_perm,
  seed = 42L
)

expect_equal(
  current = score_perm$pval,
  target = score_perm_2$pval,
  info = "singscore single - permutation reproducible with same seed"
)

### missing genes -------------------------------------------------------------

# Mix of present and missing — should still score on the present ones
mixed_set <- c(up_set[1:10], bad_set)
score_mixed <- calc_singscore(ranks = ranks, up_set = mixed_set)

expect_equal(
  current = nrow(score_mixed),
  target = no_samples,
  info = "singscore single - handles partial gene-set overlap"
)

expect_error(
  current = calc_singscore(ranks = ranks, up_set = bad_set),
  info = "singscore single - errors when no genes overlap"
)

## tests multiple gene sets ----------------------------------------------------

### up only --------------------------------------------------------------------

multi_up <- calc_singscore_multi(
  ranks = ranks,
  up_pathways = up_pathways
)

expect_true(
  current = all(c("scores", "dispersions") %in% names(multi_up)),
  info = "singscore multi - returns scores and dispersions"
)

expect_equal(
  current = dim(multi_up$scores),
  target = c(no_gs, no_samples),
  info = "singscore multi - scores shape (gene_sets x samples)"
)

expect_equal(
  current = rownames(multi_up$scores),
  target = names(up_pathways),
  info = "singscore multi - rownames = gene set names"
)

expect_equal(
  current = colnames(multi_up$scores),
  target = colnames(data),
  info = "singscore multi - colnames = sample names"
)

### up + down ------------------------------------------------------------------

multi_updn <- calc_singscore_multi(
  ranks = ranks,
  up_pathways = up_pathways,
  down_pathways = down_pathways
)

expect_equal(
  current = dim(multi_updn$scores),
  target = c(no_gs, no_samples),
  info = "singscore multi up+down - scores shape"
)

# Mismatched names: should error or drop unmatched
expect_error(
  current = calc_singscore_multi(
    ranks = ranks,
    up_pathways = up_pathways,
    down_pathways = setNames(down_pathways, sprintf("other_%i", 1:no_gs))
  ),
  info = "singscore multi - errors on no matching names"
)

## comparison against singscore ------------------------------------------------

if (!requireNamespace("singscore", quietly = TRUE)) {
  exit_file("singscore not available")
}

# Note: singscore::rankGenes uses ties.method = "min", bixverse uses average.
# With continuous Gaussian data ties are effectively absent, so results match
# tightly.

sing_ranks <- singscore::rankGenes(data)

### ranking --------------------------------------------------------------------

expect_equivalent(
  current = ranks_rs,
  target = sing_ranks,
  info = "singscore - both methods return the same ranking"
)

### single ---------------------------------------------------------------------

sing_up <- singscore::simpleScore(rankData = sing_ranks, upSet = up_set)
our_up <- calc_singscore(ranks = ranks, up_set = up_set)

expect_equivalent(
  current = our_up$total_score,
  target = sing_up$TotalScore,
  info = "singscore - up-only total score vs singscore::simpleScore",
  tolerance = 1e-9
)

expect_equivalent(
  current = our_up$total_dispersion,
  target = sing_up$TotalDispersion,
  info = "singscore - up-only dispersion vs singscore::simpleScore",
  # R's mad() implementation uses a less precise number here... Rust is MORE
  # precise
  tolerance = 1e-3
)

sing_updn <- singscore::simpleScore(
  rankData = sing_ranks,
  upSet = up_set,
  downSet = down_set
)
our_updn <- calc_singscore(
  ranks = ranks,
  up_set = up_set,
  down_set = down_set
)

expect_equivalent(
  current = our_updn$total_score,
  target = sing_updn$TotalScore,
  info = "singscore - up+down total score vs singscore::simpleScore",
  tolerance = 1e-9
)

expect_equivalent(
  current = our_updn$up_score,
  target = sing_updn$UpScore,
  info = "singscore - up component vs singscore::simpleScore",
  tolerance = 1e-9
)

expect_equivalent(
  current = our_updn$down_score,
  target = sing_updn$DownScore,
  info = "singscore - down component vs singscore::simpleScore",
  tolerance = 1e-9
)

### multi ---------------------------------------------------------------------

if (!requireNamespace("GSEABase", quietly = TRUE)) {
  exit_file("GSEABase not available")
}

up_gsc <- GSEABase::GeneSetCollection(
  lapply(names(up_pathways), function(n) {
    GSEABase::GeneSet(up_pathways[[n]], setName = n)
  })
)

sing_multi <- singscore::multiScore(rankData = sing_ranks, upSetColc = up_gsc)
our_multi <- calc_singscore_multi(
  ranks = ranks,
  up_pathways = up_pathways,
  min_size = 1L
)

common_gs <- intersect(rownames(our_multi$scores), rownames(sing_multi$Scores))

expect_equal(
  current = length(common_gs),
  target = no_gs,
  info = "singscore multi - all gene sets retained"
)

expect_equivalent(
  current = our_multi$scores[common_gs, ],
  target = sing_multi$Scores[common_gs, ],
  info = "singscore multi - scores vs singscore::multiScore",
  tolerance = 1e-9
)

expect_equivalent(
  current = our_multi$dispersions[common_gs, ],
  target = sing_multi$Dispersions[common_gs, ],
  info = "singscore multi - dispersions vs singscore::multiScore",
  # R's mad() implementation uses a less precise number here... Rust is MORE
  # precise
  tolerance = 1e-3
)
