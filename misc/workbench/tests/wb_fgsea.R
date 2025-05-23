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

# Elimination method for the gene ontology -------------------------------------

devtools::load_all()

rextendr::document()

go_data_dt <- get_go_human_data()

go_data_s7 <- gene_ontology_data(go_data_dt, min_genes = 3L)

protein_coding_genes <- unique(unlist(go_data_s7@go_to_genes))

stats <- setNames(
  sort(rnorm(length(protein_coding_genes)), decreasing = TRUE),
  protein_coding_genes
)

levels <- names(S7::prop(go_data_s7, "levels"))

tictoc::tic()
test_1 <- rs_geom_elim_fgsea_simple(
  stats = stats,
  levels = levels,
  go_obj = go_data_s7,
  gsea_param = 1.0,
  elim_threshold = 0.05,
  min_size = 5,
  max_size = 500,
  iters = 2000,
  seed = 10101,
  debug = FALSE
)
tictoc::toc()

leading_edges <- mapply(
  "[",
  list(names(stats)),
  test_1$leading_edge,
  SIMPLIFY = FALSE
)

names(test_1)


res_dt <- data.table::setDT(test_1[c(
  "go_id",
  "es",
  "nes",
  "size",
  "pvals"
)]) %>%
  .[, leading_edge := leading_edges]

test_1
