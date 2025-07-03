library(magrittr)
library(data.table)

devtools::document()
rextendr::document()
# devtools::load_all()
# devtools::check()

# Community detections ----

edge_data <- arrow::read_parquet("~/Desktop/string_clean.parquet") %>%
  as.data.table()
#
# edge_data_clean = edge_data %>%
#   .[interactionScore >= 0.85] %>%
#   dplyr::select(from = `:START_ID`, to = `:END_ID`)

devtools::load_all()

test_class <- network_diffusions(edge_data, weighted = FALSE, directed = FALSE)


set.seed(123)
genes <- sample(igraph::V(test_class@graph)$name, 10)
genes.2 <- c(sample(igraph::V(test_class@graph)$name, 10), genes[1:3])
genes.3 <- sample(igraph::V(test_class@graph)$name, 25)
diffusion_vector <- rep(1, 10) %>% `names<-`(genes)
diffusion_vector.2 <- rep(1, 13) %>% `names<-`(genes.2)

test_class <- diffuse_seed_nodes(test_class, diffusion_vector, "max")

# rust speed tests

object = test_class
perm_iters = 1000L
random_seed = 10101L
.verbose = TRUE

graph <- S7::prop(object, "graph")
diffusion_params <- S7::prop(object, "params")[["diffusion_params"]]
diffusion_results <- S7::prop(object, "diffusion_res")
nodes_names <- igraph::V(graph)$name

diffusion_vector <- diffusion_params[["diffusion_vector"]]

randomised_diffusions <- generate_perm_diffusion_vecs(
  graph = graph,
  diffusion_vec = diffusion_vector,
  iters = perm_iters
)

# prepare data for rust
edge_list <- igraph::as_edgelist(graph, names = TRUE)
graph_names <- igraph::V(graph)$name

# original

tictoc::tic()
page_rank_perm_res <- rs_page_rank_permutations(
  node_names = graph_names,
  from = edge_list[, 1],
  to = edge_list[, 2],
  diffusion_scores = randomised_diffusions,
  undirected = !igraph::is_directed(graph)
)
tictoc::toc()


# new version with claude
tictoc::tic()
page_rank_perm_res <- rs_page_rank_permutations_optimized(
  node_names = graph_names,
  from = edge_list[, 1],
  to = edge_list[, 2],
  diffusion_scores = randomised_diffusions,
  undirected = !igraph::is_directed(graph)
)
tictoc::toc()

plot(
  page_rank_perm_res$means,
  page_rank_perm_res$means
)


# test_class <- permute_seed_nodes(test_class)

test_class <- tied_diffusion(
  object = test_class,
  diffusion_vector_1 = diffusion_vector,
  diffusion_vector_2 = diffusion_vector.2,
  summarisation = "max",
  score_aggregation = "min"
)

test_class <- permute_seed_nodes(test_class)

get_params(test_class, TRUE, TRUE)

test_class <- community_detection(
  test_class,
  community_params = params_community_detection(
    min_seed_nodes = 0L,
    threshold_type = "pval_based"
  )
)

x <- get_results(test_class)

print(x)

# RBH graph --------------------------------------------------------------------

set.seed(123)

modules_set_a <- purrr::map(
  1:5,
  ~ {
    sample(letters, 10)
  }
)
names(modules_set_a) <- LETTERS[1:5]
data_a <- data.table::setDT(stack(modules_set_a))[, `:=`(
  origin = "set_A",
  ind = as.character(ind)
)]

modules_set_b <- purrr::map(
  1:3,
  ~ {
    sample(letters, 8)
  }
)
names(modules_set_b) <- LETTERS[1:3]
data_b <- data.table::setDT(stack(modules_set_b))[, `:=`(
  origin = "set_B",
  ind = as.character(ind)
)]

full_data <- rbindlist(list(data_a, data_b))

rbh_class <- rbh_graph(
  full_data,
  dataset_col = "origin",
  module_col = "ind",
  value_col = "values"
)

module_data <- rbh_class@module_data

rextendr::document()

test <- rs_rbh_sets(
  module_list = module_data,
  overlap_coefficient = TRUE,
  min_similarity = 0,
  debug = TRUE
)

rbh_class <- generate_rbh_graph(
  rbh_class,
  minimum_similarity = 0,
  overlap_coefficient = T,
  .debug = FALSE
)

rbh_class <- find_rbh_communities(rbh_class)

plot_resolution_res(rbh_class)
