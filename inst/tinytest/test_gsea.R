# (fast) gsea tests ------------------------------------------------------------

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

### gene set indices -----------------------------------------------------------

expected_pos_idx <- c(
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
expected_neg_idx <- c(855, 860, 919, 964, 972, 979, 980)
expected_rnd1_idx <- c(203, 225, 255, 354, 457, 554, 561, 577, 651, 946)
expected_rnd2_idx <- c(
  37,
  113,
  138,
  223,
  232,
  242,
  318,
  473,
  511,
  526,
  605,
  856,
  873,
  974,
  980
)
expected_rnd3_idx <- c(
  49,
  564,
  662,
  684,
  767,
  779,
  967,
  992
)

rs_indices <- rs_get_gs_indices(
  gene_universe = gene_universe,
  pathway_list = pathway_list
)

expect_equal(
  current = rs_indices$pathway_pos,
  target = expected_pos_idx,
  info = paste(
    "fgsea: positive gene set indices"
  )
)

expect_equal(
  current = rs_indices$pathway_neg,
  target = expected_neg_idx,
  info = paste(
    "fgsea: negative gene set indices"
  )
)

expect_equal(
  current = rs_indices$pathway_pos,
  target = expected_pos_idx,
  info = paste(
    "fgsea: positive gene set indices"
  )
)

expect_equal(
  current = rs_indices$random_p1,
  target = expected_rnd1_idx,
  info = paste(
    "fgsea: random gene set indices (set1)"
  )
)

expect_equal(
  current = rs_indices$random_p2,
  target = expected_rnd2_idx,
  info = paste(
    "fgsea: random gene set indices (set2)"
  )
)

expect_equal(
  current = rs_indices$random_p3,
  target = expected_rnd3_idx,
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
  current = rs_gsea_stats_pos$gene_stat,
  target = expected_gsea_stats_pos,
  info = paste(
    "fgsea: gsea stat pos"
  ),
  tolerance = 1e-6
)

expect_equal(
  current = rs_gsea_stats_pos_v2$gene_stat,
  target = expected_gsea_stats_pos,
  info = paste(
    "fgsea: gsea stat pos (with leading edge)"
  ),
  tolerance = 1e-6
)

expect_equal(
  current = rs_gsea_stats_pos_v2$leading_edge,
  target = expected_gsea_pos_leading_edge,
  info = paste(
    "fgsea: gsea stat pos: leading edge"
  ),
  tolerance = 1e-6
)

#### negative ------------------------------------------------------------------

expect_equal(
  current = rs_gsea_stats_neg$gene_stat,
  target = expected_gsea_stats_neg,
  info = paste(
    "fgsea: gsea stat neg"
  ),
  tolerance = 1e-6
)

expect_equal(
  current = rs_gsea_stats_neg$leading_edge,
  target = expected_gsea_neg_leading_edge,
  info = paste(
    "fgsea: gsea stat pos: leading edge"
  ),
  tolerance = 1e-6
)

## direct comparison fgsea vs internak -----------------------------------------

fgsea_scores <- fgsea::fgseaSimple(
  pathways = pathway_list,
  stats = stats,
  nperm = 100
)

bixverse_scores <- fgsea_simple(
  stats = stats,
  pathways = pathway_list,
  nperm = 100L
)

correlation_fgsea_internal_pval <- cor(
  x = fgsea_scores$pval,
  y = bixverse_scores$pval,
  method = "pearson"
)

correlation_fgsea_internal_es <- cor(
  x = fgsea_scores$ES,
  y = bixverse_scores$ES,
  method = "pearson"
)

correlation_fgsea_internal_nes <- cor(
  x = fgsea_scores$NES,
  y = bixverse_scores$NES,
  method = "pearson"
)

# There should be a very high correlation, despite random initialisation
expect_true(
  correlation_fgsea_internal_pval >= 0.97,
  info = paste(
    "correlation internal fgsea vs official (pval)"
  )
)

# There should be a very high correlation, despite random initialisation
expect_true(
  correlation_fgsea_internal_es >= 0.97,
  info = paste(
    "correlation internal fgsea vs official (ES)"
  )
)

# There should be a very high correlation, despite random initialisation
expect_true(
  correlation_fgsea_internal_nes >= 0.97,
  info = paste(
    "correlation internal fgsea vs official (NES)"
  )
)
