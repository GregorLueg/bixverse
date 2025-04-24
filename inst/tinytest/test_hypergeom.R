# gse tests --------------------------------------------------------------------

## simple version --------------------------------------------------------------

target_genes <- c("a", "b", "c")
gene_set <- list(gene_set_a = c("b", "c", "d"))
gene_universe <- letters

trials <- length(target_genes)
hits <- length(intersect(target_genes, gene_set[[1]]))
gene_set_length <- length(gene_set[[1]])
gene_universe_length <- length(gene_universe)

r_pval <- phyper(
  q = hits - 1,
  m = gene_set_length,
  n = gene_universe_length - gene_set_length,
  k = trials,
  lower.tail = FALSE
)

expected_result <- list(
  pvals = r_pval,
  odds_ratios = 44,
  hits = 2,
  gene_set_lengths = 3
)

rs_res <- rs_hypergeom_test(
  target_genes = target_genes,
  gene_sets = gene_set,
  gene_universe = gene_universe
)

expect_equal(
  current = rs_res,
  target = expected_result,
  info = "Gene set enrichment test with hypergeometric test"
)
