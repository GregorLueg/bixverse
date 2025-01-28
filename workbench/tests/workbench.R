library(magrittr)
library(data.table)

devtools::document()
rextendr::document()
devtools::load_all()
devtools::check()

edge_data = arrow::read_parquet(
  "~/Desktop/biomind_downloads/processed_data/edges_OT_interactions.parquet"
) %>%
  as.data.table() %>%
  .[network_resource == 'SIGNOR']

edge_data_clean = edge_data %>%
  dplyr::select(from = `:START_ID`, to = `:END_ID`)

test_class = network_diffusions(edge_data_clean, weighted = FALSE, directed = TRUE)

tictoc::tic()
calculate_diff_AUC(test_class, genes.3)
tictoc::toc()

?igraph::page_rank

set.seed(123)
genes = sample(igraph::V(test_class@graph)$name, 10)
genes.2 = sample(igraph::V(test_class@graph)$name, 10)
genes.3 = sample(igraph::V(test_class@graph)$name, 25)
diffusion_vector = rep(1, 10) %>% `names<-`(genes)
diffusion_vector.2 = rep(1, 10) %>% `names<-`(genes.2)

test_class <- diffuse_seed_genes(test_class, diffusion_vector, 'max')

calculate_diff_AUC(test_class, hit_nodes = genes.3)

test_class <- tied_diffusion(
  network_diffusions = test_class,
  diffusion_vector.1 = diffusion_vector,
  diffusion_vector.2 = diffusion_vector.2,
  summarisation = 'max',
  score_aggregation = 'min'
)

scores <- test_class@diffusion_res
pos.scores <- scores[genes.3]
neg.scores <- scores[which(!names(scores) %in% genes.3)]

?rs_fast_auc

rs_fast_auc(pos.scores, neg.scores, 10000, 42)

?rs_create_random_aucs

tictoc::tic()
random_aucs = rs_create_random_aucs(
  scores,
  size_pos = length(pos.scores),
  random_iters = 1000,
  auc_iters = 10000,
  seed = 42
)
tictoc::toc()

hist(random_aucs)
