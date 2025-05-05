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
