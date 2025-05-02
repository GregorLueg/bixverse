# Rust implementation of fgsea ----

stats <- setNames(sort(rnorm(1000), decreasing = TRUE), paste0("gene", 1:1000))
pathways <- list(
  pathway1 = paste0("gene", sample(1:100, 50)),
  pathway2 = paste0("gene", sample(1:250, 40)),
  pathway3 = paste0("gene", sample(1:1000, 100)),
  pathway5 = paste0("gene", sample(750:1000, 70)),
  pathway5 = paste0("gene", sample(900:1000, 80))
)

tictoc::tic()
results <- fgseaSimple(pathways, stats, nPermutations = 2000)
tictoc::toc()

rextendr::document()

plot(
  results$ES,
  rs_results$es
)

plot(
  results$NES,
  rs_results$nes
)

results$NES
rs_results$nes

results$pval
rs_results$pvals

rextendr::document()

devtools::load_all()

set.seed(123)

gene_size <- 20000
gene_set_size_range <- 50:100
no_gene_sets <- 500

stats <- setNames(
  sort(rnorm(gene_size), decreasing = TRUE),
  paste0("gene", 1:gene_size)
)

gene_universe <- names(stats)

random_gene_sets <- purrr::map(
  1:no_gene_sets,
  ~ {
    sample(gene_universe, sample(gene_set_size_range, 1), )
  }
)

names(random_gene_sets) <- paste0("pathway", 1:no_gene_sets)

table(purrr::map_dbl(random_gene_sets, length)) %>% sort(decreasing = TRUE)

tictoc::tic()
rs_results <- rs_gsea_traditional(
  ranks = stats,
  vec_name = names(stats),
  iters = 2000,
  pathway_list = random_gene_sets,
  seed = 123L
)
tictoc::toc()


tictoc::tic()
rs_results_2 <- rs_gsea_fgsea_simple(
  ranks = stats,
  vec_name = names(stats),
  iters = 2000,
  pathway_list = random_gene_sets,
  gsea_param = 1,
  seed = 123L
)
tictoc::toc()

plot(
  rs_results$pvals,
  rs_results_2$pvals
)

cor(rs_results$pvals, rs_results_2$pvals, method = 'spearman')

# Dissect the original implementation and identify differences ... ----

set.seed(123)

gene_size <- 20000
gene_set_size_range <- 50:100
no_gene_sets <- 500

stats <- setNames(
  sort(rnorm(gene_size), decreasing = TRUE),
  paste0("gene", 1:gene_size)
)

gene_universe <- names(stats)

random_gene_sets <- purrr::map(
  1:no_gene_sets,
  ~ {
    sample(gene_universe, sample(gene_set_size_range, 1), )
  }
)

names(random_gene_sets) <- paste0("pathway", 1:no_gene_sets)

rextendr::document()

devtools::load_all()

tictoc::tic()
internal_implementation <- fgsea_simple(
  stats = stats,
  pathways = random_gene_sets
)
tictoc::toc()


tictoc::tic()
fgsea_implementaiton <- fgsea::fgseaSimple(
  pathways = random_gene_sets,
  stats = stats,
  nperm = 2000
)
tictoc::toc()

plot(internal_implementation$pval, fgsea_implementaiton$pval)
