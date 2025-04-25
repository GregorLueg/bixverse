# gse tests --------------------------------------------------------------------

## simple version --------------------------------------------------------------

### data -----------------------------------------------------------------------

target_genes <- c("a", "b", "c")
target_genes_list <- list(
  first_set = target_genes,
  second_set = c("d", "e")
)
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

expected_results_dt <- data.table::as.data.table(
  expected_result
) %>%
  .[, `:=`(
    gene_set_name = "gene_set_a",
    fdr = pvals,
    target_set_lengths = trials
  )] %>%
  .[,
    c(
      "gene_set_name",
      "odds_ratios",
      "pvals",
      "fdr",
      "hits",
      "gene_set_lengths",
      "target_set_lengths"
    ),
    with = FALSE
  ]

### rust vs rs -----------------------------------------------------------------

rs_res <- rs_hypergeom_test(
  target_genes = target_genes,
  gene_sets = gene_set,
  gene_universe = gene_universe
)

expect_equal(
  current = rs_res,
  target = expected_result,
  info = paste(
    "Gene set enrichment test with hypergeometric test:",
    "(Rust <> R equivalence)"
  )
)

### direct function ------------------------------------------------------------

#### single --------------------------------------------------------------------

r_res <- gse_hypergeometric(
  target_genes = target_genes,
  gene_set_list = gene_set,
  gene_universe = gene_universe,
  threshold = 1,
  minimum_overlap = 1L
)

expect_equal(
  current = r_res,
  target = expected_results_dt,
  info = paste(
    "Gene set enrichment test with hypergeometric test:",
    "(Wrapper function)"
  )
)

#### multiple ------------------------------------------------------------------

expected_res_multiple <- data.table::data.table(
  target_set_name = names(target_genes_list),
  odds_ratios = c(44, 11),
  pvals = c(r_pval, 1),
  fdr = c(r_pval, 1),
  hits = c(2, 1),
  gene_set_lengths = 3,
  gene_set_name = names(gene_set),
  target_set_lengths = c(3, 2)
)

res_multiple <- gse_hypergeometric_list(
  target_genes_list = target_genes_list,
  gene_set_list = gene_set,
  gene_universe = gene_universe,
  threshold = 1,
  minimum_overlap = 1L
)

expect_equal(
  current = expected_res_multiple,
  target = res_multiple,
  info = paste(
    "Gene set enrichment test with hypergeometric test:",
    "(multiple tests)"
  )
)

## gene ontology elimination methods -------------------------------------------

### data -----------------------------------------------------------------------
