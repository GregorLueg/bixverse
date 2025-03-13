library(magrittr)
library(data.table)

devtools::document()
rextendr::document()
devtools::load_all()
devtools::check()

# Community detections ----

edge_data = arrow::read_parquet("~/Desktop/biomind_downloads/processed_data/edges_interactions.parquet") %>%
  as.data.table() %>%
  .[network_resource == 'STRING']

edge_data_clean = edge_data %>%
  .[interactionScore >= 0.85] %>%
  dplyr::select(from = `:START_ID`, to = `:END_ID`)

test_class = network_diffusions(edge_data_clean, weighted = FALSE, directed = TRUE)


set.seed(123)
genes = sample(igraph::V(test_class@graph)$name, 10)
genes.2 = sample(igraph::V(test_class@graph)$name, 10)
genes.3 = sample(igraph::V(test_class@graph)$name, 25)
diffusion_vector = rep(1, 10) %>% `names<-`(genes)
diffusion_vector.2 = rep(1, 10) %>% `names<-`(genes.2)


test_class <- diffuse_seed_nodes(test_class, diffusion_vector, 'max')

get_params(test_class, TRUE, TRUE)

devtools::load_all()

test_class <- tied_diffusion(
  object = test_class,
  diffusion_vector_1 = diffusion_vector,
  diffusion_vector_2 = diffusion_vector.2,
  summarisation = 'max',
  score_aggregation = 'min'
)


get_params(test_class, TRUE, TRUE)

test_class <- community_detection(
  test_class,
  min_seed_nodes = 0L,
  diffusion_threshold = .5,
  max_nodes = 300L
)

get_params(test_class, to_json = TRUE, pretty_json = TRUE)

x = get_results(test_class)

?get_results

test_class@params

devtools::load_all()

get_params(test_class, TRUE, TRUE)

?jsonlite::toJSON

jsonlite::toJSON(test_class@params) %>% jsonlite::prettify()

# RBH graph ----

protein_coding_genes <- data.table::fread("~/Desktop/protein_coding_genes.csv")

universe <- protein_coding_genes$id[1:500]

sets_per_origin <- 100
gene_sets_no <- 100

gene_sets_no <- sets_per_origin * length(LETTERS)

seed <- 123
set.seed(seed)

gene_sets <- purrr::map(1:gene_sets_no, ~ {
  set.seed(seed + .x + 1)
  size <- sample(20:100, 1)
  sample(universe, size, replace = FALSE)
})

names(gene_sets) <- purrr::map_chr(1:gene_sets_no, ~ {
  set.seed(seed + .x + 1)
  paste(sample(LETTERS, 5), collapse = "")
})

gene_sets_df <- purrr::imap(gene_sets, ~{
  data.table::data.table(
    genes = .x,
    name = .y
  )
})

origins <- rep(LETTERS, each = sets_per_origin)

gene_sets_df <- purrr::map2(gene_sets_df, origins, ~{
  .x[, origin := .y]
}) %>% rbindlist


module_df = gene_sets_df

devtools::document()

rbh_class = rbh_graph(
  gene_sets_df,
  dataset_col = 'origin',
  module_col = 'name',
  value_col = 'genes'
)


rbh_class = generate_rbh_graph(rbh_class, minimum_similarity = 0.2, overlap_coefficient = T, .debug = FALSE)

object <- rbh_class

list_of_list <- S7::prop(object, "module_data")

list_of_list <- list(
  "dataset_A" = list(
    "module_A" = c("A", "B", "C"),
    "module_B" = c("D", "E", "H")
  ),
  "dataset_B" = list(
    "module_A" = c("A", "B", "C"),
    "module_B" = c("X", "A", "Z")
  ),
  "dateset_B" = list(
    "test_A" = c("A", "B")
  )
)

rbh_results <- rs_rbh_sets(
  module_list = list_of_list,
  overlap_coefficient = FALSE,
  min_similarity = 0.2,
  debug = TRUE
)




rextendr::document()

rbh_results$similarity

rbh_class@rbh_edge_df

rbh_class = find_rbh_communities(rbh_class)

get_params(rbh_class, TRUE, TRUE)

get_results(rbh_class)

rbh_class

edge_dt = S7::prop(rbh_class, "rbh_edge_df")

graph = S7::prop(rbh_class, "rbh_graph")



cluster_red <- igraph::cluster_louvain(
  graph,
  resolution = .5
)

?igraph::cluster_infomap

clusters_red <- igraph::cluster_infomap(
  graph
)

community_res <- data.table(
  node_id = igraph::V(graph)$name,
  community_name = clusters_red$membership
)





table(clusters_red$membership)

graph
