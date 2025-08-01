# (fast) gsea tests ------------------------------------------------------------

library(magrittr)

set.seed(123)

stat_size <- 1000

stats <- setNames(
  sort(rnorm(stat_size), decreasing = TRUE),
  paste0("gene", 1:stat_size)
)

random_sizes <- c(10, 15, 8)
pathway_random <- purrr::map(
  random_sizes,
  ~ {
    sample(names(stats), .x)
  }
)
pathway_pos <- sample(names(stats)[1:150], 15)
pathway_neg <- sample(names(stats)[851:1000], 7)
gene_universe <- names(stats)

pathway_list <- list(
  pathway_pos = pathway_pos,
  pathway_neg = pathway_neg,
  random_p1 = pathway_random[[1]],
  random_p2 = pathway_random[[2]],
  random_p3 = pathway_random[[3]]
)

## general functionality -------------------------------------------------------

### enrichment scores ----------------------------------------------------------

rs_es_pos <- rs_calc_es(
  stats = stats,
  pathway_r = pathway_pos
)

rs_es_neg <- rs_calc_es(
  stats = stats,
  pathway_r = pathway_neg
)

expect_equal(
  current = rs_es_pos,
  target = 0.8730964,
  info = paste(
    "gsea: es calculation positive"
  ),
  tolerance = 10e-7
)

expect_equal(
  current = rs_es_neg,
  target = -0.8600201,
  info = paste(
    "gsea: es calculation negative"
  ),
  tolerance = 10e-7
)

### gene set indices -----------------------------------------------------------

pathway_indices_r <- purrr::map(pathway_list, \(genes) {
  which(gene_universe %in% genes)
})

rs_indices <- rs_get_gs_indices(
  gene_universe = gene_universe,
  pathway_list = pathway_list
)

expect_equal(
  current = rs_indices$pathway_pos,
  target = pathway_indices_r$pathway_pos,
  info = paste(
    "fgsea: positive gene set indices"
  )
)

expect_equal(
  current = rs_indices$pathway_neg,
  target = pathway_indices_r$pathway_neg,
  info = paste(
    "fgsea: negative gene set indices"
  )
)

expect_equal(
  current = rs_indices$random_p1,
  target = pathway_indices_r$random_p1,
  info = paste(
    "fgsea: random gene set indices (set1)"
  )
)

expect_equal(
  current = rs_indices$random_p2,
  target = pathway_indices_r$random_p2,
  info = paste(
    "fgsea: random gene set indices (set2)"
  )
)

expect_equal(
  current = rs_indices$random_p3,
  target = pathway_indices_r$random_p3,
  info = paste(
    "fgsea: random gene set indices (set3)"
  )
)

### gene set statistics --------------------------------------------------------

expected_gsea_stats_pos <- 0.8730964

expected_gsea_stats_neg <- -0.8600201

expected_gsea_pos_leading_edge <- c(
  18,
  30,
  32,
  36,
  37,
  46,
  49,
  51,
  70,
  82,
  110,
  112,
  127,
  131,
  140
)
expected_gsea_neg_leading_edge <- c(
  980,
  979,
  972,
  964,
  919,
  860,
  855
)

rs_gsea_stats_pos = rs_calc_gsea_stats(
  stats = stats,
  gs_idx = rs_indices$pathway_pos,
  gsea_param = 1.0,
  return_leading_edge = FALSE
)

rs_gsea_stats_pos_v2 = rs_calc_gsea_stats(
  stats = stats,
  gs_idx = rs_indices$pathway_pos,
  gsea_param = 1.0,
  return_leading_edge = TRUE
)

rs_gsea_stats_neg = rs_calc_gsea_stats(
  stats = stats,
  gs_idx = rs_indices$pathway_neg,
  gsea_param = 1.0,
  return_leading_edge = TRUE
)


#### positive ------------------------------------------------------------------

expect_equal(
  current = rs_gsea_stats_pos$leading_edge,
  target = vector(mode = "integer"),
  info = paste(
    "gsea: no leading edge genes returned"
  )
)

expect_equal(
  current = rs_gsea_stats_pos$es,
  target = expected_gsea_stats_pos,
  info = paste(
    "gsea: gsea stat pos"
  ),
  tolerance = 1e-6
)

expect_equal(
  current = rs_gsea_stats_pos_v2$es,
  target = expected_gsea_stats_pos,
  info = paste(
    "gsea: gsea stat pos (with leading edge)"
  ),
  tolerance = 1e-7
)

expect_equal(
  current = rs_gsea_stats_pos_v2$leading_edge,
  target = expected_gsea_pos_leading_edge,
  info = paste(
    "gsea: gsea stat pos: leading edge"
  ),
  tolerance = 1e-7
)

#### negative ------------------------------------------------------------------

expect_equal(
  current = rs_gsea_stats_neg$es,
  target = expected_gsea_stats_neg,
  info = paste(
    "gsea: gsea stat neg"
  ),
  tolerance = 1e-7
)

expect_equal(
  current = rs_gsea_stats_neg$leading_edge,
  target = expected_gsea_neg_leading_edge,
  info = paste(
    "gsea: gsea stat pos: leading edge"
  ),
  tolerance = 1e-7
)

## traditional gsea vs simple fgsea --------------------------------------------

# generally speaking this should yield the same

internal_gsea_simple_res <- calc_fgsea_simple(
  stats = stats,
  pathways = pathway_list,
  nperm = 100L
) %>%
  data.table::setorder(pathway_name)

traditional_gsea_results <- calc_gsea_traditional(
  stats = stats,
  pathways = pathway_list,
  nperm = 100L
) %>%
  data.table::setorder(pathway_name)

correlation_traditional_vs_fgsea_es <- cor(
  internal_gsea_simple_res$es,
  traditional_gsea_results$es
)

correlation_traditional_vs_fgsea_nes <- cor(
  internal_gsea_simple_res$nes,
  traditional_gsea_results$nes
)

correlation_traditional_vs_fgsea_pval <- cor(
  internal_gsea_simple_res$pval,
  traditional_gsea_results$pvals
)

expect_true(
  correlation_traditional_vs_fgsea_es >= 0.97,
  info = paste(
    "correlation internal fgsea vs official (ES)"
  )
)

expect_true(
  correlation_traditional_vs_fgsea_nes >= 0.97,
  info = paste(
    "correlation internal fgsea vs official (NES)"
  )
)

expect_true(
  correlation_traditional_vs_fgsea_pval >= 0.97,
  info = paste(
    "correlation internal fgsea vs official (pval)"
  )
)

## multi-level version ---------------------------------------------------------

### comparison 1 ---------------------------------------------------------------

# For the two significant pathways, we should see lower p-values

internal_gsea_multi_level_res <- calc_fgsea(
  stats = stats,
  pathways = pathway_list,
  nperm = 100L
) %>%
  data.table::setorder(pathway_name)

expect_true(
  all(
    internal_gsea_multi_level_res$pvals[1:2] <
      internal_gsea_simple_res$pvals[1:2]
  ),
  info = paste(
    "simple fgsea vs multi-level fgsea comparison (low permutations)"
  )
)

### comparison 2 ---------------------------------------------------------------

# with higher number of permutations, the p-values for the random pathways
# should be the same

internal_gsea_simple_res_v2 <- calc_fgsea_simple(
  stats = stats,
  pathways = pathway_list,
  nperm = 1000L
) %>%
  data.table::setorder(pathway_name)

internal_gsea_multi_level_res_v2 <- calc_fgsea(
  stats = stats,
  pathways = pathway_list,
  nperm = 1000L
) %>%
  data.table::setorder(pathway_name)

matching_pvalues <- sum(
  internal_gsea_simple_res_v2$pvals == internal_gsea_multi_level_res_v2$pvals
)

expect_true(
  matching_pvalues == 3,
  info = paste(
    "simple fgsea vs multi-level fgsea comparison (high permutations)"
  )
)

# direct comparison fgsea vs internal ------------------------------------------

## calc gsea stats -------------------------------------------------------------

if (requireNamespace("fgsea", quietly = TRUE)) {
  fgsea_result_calc_gsea_stats_pos <- fgsea::calcGseaStat(
    stats = stats,
    selectedStats = pathway_indices_r$pathway_pos,
    returnLeadingEdge = TRUE
  ) %>%
    `names<-`(c("es", "leading_edge"))

  fgsea_result_calc_gsea_stats_neg <- fgsea::calcGseaStat(
    stats = stats,
    selectedStats = pathway_indices_r$pathway_neg,
    returnLeadingEdge = TRUE
  ) %>%
    `names<-`(c("es", "leading_edge"))

  expect_equal(
    current = rs_gsea_stats_pos_v2,
    target = fgsea_result_calc_gsea_stats_pos,
    info = paste(
      "calc gsea stats fgsea vs internal: positive"
    )
  )

  expect_equal(
    current = rs_gsea_stats_neg,
    target = fgsea_result_calc_gsea_stats_neg,
    info = paste(
      "calc gsea stats fgsea vs internal: positive"
    )
  )
}

## simple method ---------------------------------------------------------------

# Check if fgsea is installed
if (requireNamespace("fgsea", quietly = TRUE)) {
  fgsea_scores <- fgsea::fgseaSimple(
    pathways = pathway_list,
    stats = stats,
    nperm = 100
  ) %>%
    data.table::setorder(pathway)

  correlation_fgsea_internal_pval <- cor(
    x = log(fgsea_scores$pval),
    y = log(internal_gsea_simple_res$pval),
    method = "pearson"
  )
  correlation_fgsea_internal_es <- cor(
    x = fgsea_scores$ES,
    y = internal_gsea_simple_res$es,
    method = "pearson"
  )
  correlation_fgsea_internal_nes <- cor(
    x = fgsea_scores$NES,
    y = internal_gsea_simple_res$nes,
    method = "pearson"
  )

  # There should be a very high correlation, despite random initialisation
  expect_true(
    correlation_fgsea_internal_pval >= 0.99,
    info = paste(
      "correlation internal fgsea vs official (pval)"
    )
  )
  # There should be a very high correlation, despite random initialisation
  expect_true(
    correlation_fgsea_internal_es >= 0.99,
    info = paste(
      "correlation internal fgsea vs official (ES)"
    )
  )
  # There should be a very high correlation, despite random initialisation
  expect_true(
    correlation_fgsea_internal_nes >= 0.99,
    info = paste(
      "correlation internal fgsea vs official (NES)"
    )
  )
}

## multi-level method ----------------------------------------------------------

if (requireNamespace("fgsea", quietly = TRUE)) {
  fgsea_scores_ml <- fgsea::fgsea(
    pathways = pathway_list,
    stats = stats
  ) %>%
    data.table::setorder(pathway)

  correlation_fgsea_multilevel_pval <- cor(
    x = log(fgsea_scores_ml$pval),
    y = log(internal_gsea_multi_level_res_v2$pvals),
    method = "pearson"
  )

  correlation_fgsea_multilevel_es <- cor(
    x = fgsea_scores_ml$ES,
    y = internal_gsea_multi_level_res_v2$es,
    method = "pearson"
  )

  correlation_fgsea_multilevel_nes <- cor(
    x = fgsea_scores_ml$NES,
    y = internal_gsea_multi_level_res_v2$nes,
    method = "pearson"
  )

  # There should be a very high correlation, despite random initialisation
  expect_true(
    correlation_fgsea_multilevel_pval >= 0.99,
    info = paste(
      "correlation internal fgsea vs official - multilevel (pval)"
    )
  )
  # There should be a very high correlation, despite random initialisation
  expect_true(
    correlation_fgsea_multilevel_es >= 0.99,
    info = paste(
      "correlation internal fgsea vs official - multilevel (ES)"
    )
  )
  # There should be a very high correlation, despite random initialisation
  expect_true(
    correlation_fgsea_multilevel_nes >= 0.99,
    info = paste(
      "correlation internal fgsea vs official - multilevel (NES)"
    )
  )
}
