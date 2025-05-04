# Rust implementation of fgsea ----

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

rextendr::document()

## Dissect the function ----

stats = stats
pathways = pathway_list
nperm = 2000L
gsea_params = params_gsea()
seed = 123L

checkmate::assertNumeric(stats, min.len = 3L, finite = TRUE)
checkmate::assertNames(names(stats))
checkmate::assertList(pathways, types = "character")
checkmate::assertNames(names(pathways))
assertGSEAParams(gsea_params)

c(stats, pathways_clean, pathway_sizes) %<-%
  with(
    gsea_params,
    prep_stats_pathways(
      stats = stats,
      pathways = pathways,
      min_size = min_size,
      max_size = max_size
    )
  )

gsea_stat_res <- with(
  gsea_params,
  do.call(
    rbind,
    lapply(
      pathways_clean,
      rs_calc_gsea_stats,
      stats = stats,
      gsea_param = gsea_param,
      return_leading_edge = TRUE
    )
  )
)

leading_edges <- mapply(
  "[",
  list(names(stats)),
  gsea_stat_res[, "leading_edge"],
  SIMPLIFY = FALSE
)

pathway_scores <- unlist(gsea_stat_res[, "es"])

permutations_res <- with(
  gsea_params,
  rs_calc_gsea_stat_cumulative_batch(
    stats = stats,
    pathway_scores = pathway_scores,
    pathway_sizes = as.integer(pathway_sizes),
    iters = nperm,
    seed = seed,
    gsea_param = gsea_param
  )
)

le_zero_mean <- permutations_res$le_zero_sum / permutations_res$le_zero
ge_zero_mean <- permutations_res$ge_zero_sum / permutations_res$ge_zero

nes <- data.table::fifelse(
  (pathway_scores > 0 & ge_zero_mean != 0) |
    (pathway_scores < 0 & le_zero_mean != 0),
  pathway_scores /
    data.table::fifelse(
      pathway_scores > 0,
      ge_zero_mean,
      abs(le_zero_mean)
    ),
  NA
)

pvals <- pmin(
  (1 + permutations_res$le_es) / (1 + permutations_res$le_zero),
  (1 + permutations_res$ge_es) / (1 + permutations_res$ge_zero)
)

final_res <- data.table::as.data.table(
  gsea_stat_res,
  keep.rownames = "pathway"
) %>%
  .[, `:=`(nes = nes, pval = pvals, leading_edge = leading_edges)]
