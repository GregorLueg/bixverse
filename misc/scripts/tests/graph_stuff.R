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
genes.2 <- sample(igraph::V(test_class@graph)$name, 10)
genes.3 <- sample(igraph::V(test_class@graph)$name, 25)
diffusion_vector <- rep(1, 10) %>% `names<-`(genes)
diffusion_vector.2 <- rep(1, 10) %>% `names<-`(genes.2)

test_class <- diffuse_seed_nodes(test_class, diffusion_vector, "max")

# write a permutation test...

object <- test_class
perm_iter <- 1000L
cache <- TRUE

graph <- S7::prop(object, "graph")
diffusion_params <- S7::prop(object, "params")[["diffusion_params"]]
diffusion_results <- S7::prop(object, "diffusion_res")

no_nodes <- length(igraph::V(graph))

nodes_names <- igraph::V(graph)$name
diffusion_vector <- diffusion_params[["diffusion_vector"]]

# create the randomised diffusion vectors accounting for node degree

node_degree_distribution <- log(igraph::degree(graph))

node_degree_discrete <- cut(node_degree_distribution, 25L) %>%
  `names<-`(names(node_degree_distribution))

diffusion_names <- names(diffusion_vector)

degree_groups <- split(names(node_degree_discrete), node_degree_discrete)
node_degrees <- node_degree_discrete[diffusion_names]

randomised_sets <- purrr::map(1:perm_iter, \(i) {
  set.seed(123L + i)

  random_set_i <- purrr::map_chr(node_degrees, \(degree) {
    sample(degree_groups[[as.character(degree)]], 1)
  })

  diffusion_vec_i <- diffusion_vector %>% `names<-`(random_set_i)

  seed_nodes_i <- intersect(names(diffusion_vec_i), nodes_names)

  diff_vec_i <- rep(0, length(nodes_names)) %>% `names<-`(nodes_names)
  for (node in seed_nodes_i) {
    diff_vec_i[node] <- diffusion_vec_i[node]
  }

  diff_vec_i <- diff_vec_i / sum(diff_vec_i)

  diff_vec_i
})

# Rust version ... ?

edge_list <- igraph::as_edgelist(graph, names = TRUE)
graph_names <- igraph::V(graph)$name

rs_test <- rs_page_rank_permutations(
  node_names = graph_names,
  from = edge_list[, 1],
  to = edge_list[, 2],
  diffusion_scores = randomised_sets,
  undirected = TRUE
)

z_scores <- (diffusion_results - rs_test$means) / (rs_test$sd + 10^-32)

pvals <- pnorm(abs(z_scores), lower.tail = F)

hist(pvals)

table(pvals <= 0.05)

table(z_scores > 0)

hist(z_scores)

hist(node_degree_distribution)

get_params(test_class, TRUE, TRUE)

devtools::load_all()

test_class <- tied_diffusion(
  object = test_class,
  diffusion_vector_1 = diffusion_vector,
  diffusion_vector_2 = diffusion_vector.2,
  summarisation = "max",
  score_aggregation = "min"
)

get_results(test_class)

get_params(test_class, TRUE, TRUE)

test_class <- community_detection(
  test_class,
  community_params = params_community_detection(
    min_seed_nodes = 1L
  ),
  diffusion_threshold = .5
)

x <- get_results(test_class)

x

# RBH graph ----

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
