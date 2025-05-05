# Rust implementation of fgsea ----

set.seed(123)

stat_size <- 20000

stats <- setNames(
  sort(rnorm(stat_size), decreasing = TRUE),
  paste0("gene", 1:stat_size)
)

number_gene_sets <- 5000
min_size <- 50
max_size <- 250

pathway_random <- purrr::map(
  seq_len(number_gene_sets),
  ~ {
    sample(names(stats), sample(min_size:max_size, 1))
  }
)

names(pathway_random) <- paste0("pathway", 1:number_gene_sets)

devtools::load_all()

tictoc::tic()
results_traditional <- calc_gsea_traditional(
  stats = stats,
  pathways = pathway_random
)
tictoc::toc()

tictoc::tic()
results_simple_fgsea <- calc_fgsea_simple(
  stats = stats,
  pathways = pathway_random
)
tictoc::toc()


tictoc::tic()
fgsea_scores_original <- fgsea::fgseaSimple(
  pathways = pathway_random,
  stats = stats,
  nperm = 2000L
)
tictoc::toc()

plot(results_traditional$es, results_simple_fgsea$es)

plot(results_traditional$pvals, results_simple_fgsea$pvals)

plot(fgsea_scores_original$pval, results_simple_fgsea$pvals)


## Debugging after life times ----

devtools::load_all()

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


traditional_gsea_results <- calc_gsea_traditional(
  stats = stats,
  pathways = pathway_list,
  nperm = 100L
)


internal_gsea_simple_res <- calc_fgsea_simple(
  stats = stats,
  pathways = pathway_list,
  nperm = 100L
)


pathways = pathway_list
nperm = 100L
gsea_params = params_gsea()
seed = 123L

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


rextendr::document()

permutations_res_simple <- with(
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
