# gse tests --------------------------------------------------------------------

library(magrittr)
. <- NULL

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

go_target_genes <- list(
  first_test = letters[2:4],
  second_test = letters[6:10]
)

toy_go_data <- data.table::data.table(
  go_id = sprintf("go_%i", 1:3),
  go_name = sprintf("go_name_%s", letters[1:3]),
  ancestors = list(
    c("go_1"),
    c("go_1", "go_2"),
    c("go_1", "go_2", "go_3")
  ),
  ensembl_id = list(c(letters[1:10]), c(letters[c(1:6, 11)]), c(letters[1:3])),
  depth = c(1, 2, 3)
) %>%
  data.table::setorder(-depth)

object <- gene_ontology_data(toy_go_data, min_genes = 1L)

#### scenario 1 data -----------------------------------------------------------

expected_pval_no_elim_v1 <- c(0.7272727, 0.2121212, 0.1515152)
expected_hits_no_elim_v1 <- c(3, 3, 2)

expected_pval_with_elim_v1 <- c(1, 1, 0.1515152)
expected_hits_with_elim_v1 <- c(1, 1, 2)

#### scenario 2 data -----------------------------------------------------------

# Should be the same, because no elimination should happen
expected_pval_no_elim_v2 <- expected_pval_with_elim_v2 <- c(0.5454545, 1, 1)
expected_hits_no_elim_v2 <- expected_hits_with_elim_v2 <- c(5, 1, 0)

### without elimination --------------------------------------------------------

#### scenario 1 ----------------------------------------------------------------

go_results_no_elim_v1 <- gse_go_elim_method(
  object = object,
  target_genes = go_target_genes$first_test,
  minimum_overlap = 0L,
  fdr_threshold = 1,
  elim_threshold = 0,
  .debug = FALSE
) %>%
  data.table::setorder(go_id)

expect_equal(
  current = go_results_no_elim_v1$hits,
  target = expected_hits_no_elim_v1,
  info = paste(
    "Gene ontology elimination GSE: hits, no elim (v1)"
  )
)

expect_equal(
  current = go_results_no_elim_v1$pvals,
  target = expected_pval_no_elim_v1,
  info = paste(
    "Gene ontology elimination GSE: pvals, no elim (v1)"
  ),
  tolerance = 10e-6
)

#### scenario 2 ----------------------------------------------------------------

go_results_no_elim_v2 <- gse_go_elim_method(
  object = object,
  target_genes = go_target_genes$second_test,
  minimum_overlap = 0L,
  fdr_threshold = 1,
  elim_threshold = 0,
  .debug = FALSE
) %>%
  data.table::setorder(go_id)


expect_equal(
  current = go_results_no_elim_v2$hits,
  target = expected_hits_no_elim_v2,
  info = paste(
    "Gene ontology elimination GSE: hits, no elim (v2)"
  )
)

expect_equal(
  current = go_results_no_elim_v2$pvals,
  target = expected_pval_no_elim_v2,
  info = paste(
    "Gene ontology elimination GSE: pvals, no elim (v2)"
  ),
  tolerance = 10e-6
)

### with elimination -----------------------------------------------------------

#### scenario 1 ----------------------------------------------------------------

go_results_with_elim_v1 <- gse_go_elim_method(
  object = object,
  target_genes = go_target_genes$first_test,
  minimum_overlap = 0L,
  fdr_threshold = 1,
  elim_threshold = 0.95,
  .debug = FALSE
) %>%
  data.table::setorder(go_id)

expect_equal(
  current = go_results_with_elim_v1$hits,
  target = expected_hits_with_elim_v1,
  info = paste(
    "Gene ontology elimination GSE: hits, with elim (v1)"
  )
)

expect_equal(
  current = go_results_with_elim_v1$pvals,
  target = expected_pval_with_elim_v1,
  info = paste(
    "Gene ontology elimination GSE: pvals, with elim (v1)"
  ),
  tolerance = 10e-6
)

#### scenario 2 ----------------------------------------------------------------

go_results_with_elim_v2 <- gse_go_elim_method(
  object = object,
  target_genes = go_target_genes$second_test,
  minimum_overlap = 0L,
  fdr_threshold = 1,
  elim_threshold = 0.95,
  .debug = FALSE
) %>%
  data.table::setorder(go_id)

expect_equal(
  current = go_results_with_elim_v2$hits,
  target = expected_hits_with_elim_v2,
  info = paste(
    "Gene ontology elimination GSE: hits, with elim (v2)"
  )
)

expect_equal(
  current = go_results_with_elim_v2$pvals,
  target = expected_pval_with_elim_v2,
  info = paste(
    "Gene ontology elimination GSE: pvals, with elim (v2)"
  ),
  tolerance = 10e-6
)

### multiple gene sets ---------------------------------------------------------

#### with elimination ----------------------------------------------------------

expected_pval_multi_elim <- c(0.1515152, 1, 1, 1, 1, 0.5454545)
expected_hits_multi_elim <- c(2, 1, 1, 0, 1, 5)

go_results_with_multiple_elim <- gse_go_elim_method_list(
  object = object,
  target_gene_list = go_target_genes,
  minimum_overlap = 0L,
  fdr_threshold = 1,
  elim_threshold = 0.95
) %>%
  data.table::setorder(target_set_name, -go_id)

expect_equal(
  current = go_results_with_multiple_elim$hits,
  target = expected_hits_multi_elim,
  info = paste(
    "Gene ontology elimination GSE - list with elim (hits)"
  )
)

expect_equal(
  current = go_results_with_multiple_elim$pvals,
  target = expected_pval_multi_elim,
  info = paste(
    "Gene ontology elimination GSE - list with elim (pvals)"
  ),
  tolerance = 10e-6
)


#### without elimination -------------------------------------------------------

expected_pval_multi_no_elim <- c(
  0.1515152,
  0.2121212,
  0.7272727,
  1,
  1,
  0.5454545
)
expected_hits_multi_no_elim <- c(2, 3, 3, 0, 1, 5)

go_results_with_multiple <- gse_go_elim_method_list(
  object = object,
  target_gene_list = go_target_genes,
  minimum_overlap = 0L,
  fdr_threshold = 1,
  elim_threshold = 0,
  .debug = FALSE
) %>%
  data.table::setorder(target_set_name, -go_id)

expect_equal(
  current = go_results_with_multiple$hits,
  target = expected_hits_multi_no_elim,
  info = paste(
    "Gene ontology elimination GSE - list no elim (hits)"
  )
)

expect_equal(
  current = go_results_with_multiple$pvals,
  target = expected_pval_multi_no_elim,
  info = paste(
    "Gene ontology elimination GSE - list no elim (pvals)"
  ),
  tolerance = 10e-6
)
