library(magrittr)
library(data.table)

devtools::document()
rextendr::document()
devtools::load_all()
devtools::check()


edge_data = arrow::read_parquet("~/Desktop/biomind_downloads/processed_data/edges_OT_interactions.parquet") %>%
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

test_class <- tied_diffusion(
  network_diffusions = test_class,
  diffusion_vector.1 = diffusion_vector,
  diffusion_vector.2 = diffusion_vector.2,
  summarisation = 'max',
  score_aggregation = 'min'
)

?tied_diffusion

?get_params

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

