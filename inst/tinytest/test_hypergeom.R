# gse tests --------------------------------------------------------------------

library(magrittr)

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

go_target_genes <- letters[2:4]

toy_go_data <- data.table::data.table(
  go_id = sprintf("go_%i", 1:3),
  go_name = sprintf("go_name_%s", letters[1:3]),
  ancestors = list(
    c("go_1"),
    c("go_1", "go_2"),
    c("go_1", "go_2", "go_3")
  ),
  ensembl_id = list(c(letters[1:10]), c(letters[1:6]), c(letters[1:3])),
  depth = c(1, 2, 3)
) %>%
  data.table::setorder(-depth)

object <- gene_ontology_data(toy_go_data, min_genes = 1L)

expected_pval_no_elim <- c(1, 0.1666667, 0.1833333)
expected_hits_no_elim <- c(3, 3, 2)

expected_pval_with_elim <- c(1, 1, 0.1833333)
expected_hits_with_elim <- c(0, 1, 2)

### without elimination --------------------------------------------------------

go_results_no_elim <- gse_go_elim_method(
  object = object,
  target_genes = go_target_genes,
  minimum_overlap = 0L,
  fdr_threshold = 1,
  elim_threshold = 0,
  .debug = FALSE
) %>%
  setorder(go_id)

expect_equal(
  current = go_results_no_elim$hits,
  target = expected_hits_no_elim,
  info = paste(
    "Gene ontology elimination GSE: no hits, no elimination"
  )
)

expect_equal(
  current = go_results_no_elim$pvals,
  target = expected_pval_no_elim,
  info = paste(
    "Gene ontology elimination GSE: pvals, no elimination"
  ),
  tolerance = 10e-6
)

### with elimination -----------------------------------------------------------

go_results_with_elim <- gse_go_elim_method(
  object = object,
  target_genes = go_target_genes,
  minimum_overlap = 0L,
  fdr_threshold = 1,
  elim_threshold = 1,
  .debug = FALSE
) %>%
  data.table::setorder(go_id)


expect_equal(
  current = go_results_with_elim$hits,
  target = expected_hits_with_elim,
  info = paste(
    "Gene ontology elimination GSE: no hits, no elimination"
  )
)

expect_equal(
  current = go_results_with_elim$pvals,
  target = expected_pval_with_elim,
  info = paste(
    "Gene ontology elimination GSE: pvals, no elimination"
  ),
  tolerance = 10e-6
)
