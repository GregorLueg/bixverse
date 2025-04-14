library(ontologyIndex)
library(magrittr)
library(data.table)

data("hpo")

devtools::load_all()

hpo_data <- as.data.table(stack(hpo$children)) %>%
  setnames(old = c('values', 'ind'), new = c('child', 'parent'))

hpo_onto <- ontology(hpo_data)

hpo_onto <- calculate_semantic_sim_onto(hpo_onto, sim_type = "combined")

resolution_params = params_graph_resolution(max_res = 20)
parallel = TRUE
max_workers = 5L
random_seed = 42L
.verbose = TRUE

edge_dt <- hpo_onto@semantic_similarities %>%
  `colnames<-`(c("from", "to", "weight")) %>%
  .[weight >= 0.7]

sims = edge_dt$weight

hist(sims)

affinity <- 1 - sims

hist(affinity)

sparsed <- rs_rbf_function(affinity, epsilon = 2, rbf_type = "bump")

hist(sparsed)

new_edge_dt = copy(edge_dt)[, weight := sparsed]

graph <- igraph::graph_from_data_frame(edge_dt, directed = FALSE)

igraph::is_weighted(graph)

if (parallel) {
  if (.verbose)
    message(sprintf("Using parallel computation over %i cores.", max_workers))

  # future plan funkiness
  assign(".temp_workers", max_workers, envir = .GlobalEnv)
  on.exit(rm(".temp_workers", envir = .GlobalEnv))

  plan(future::multisession(workers = .temp_workers))
} else {
  if (.verbose)
    message("Using sequential computation.")
  future::plan(future::sequential())
}

resolutions <- with(resolution_params, exp(seq(log(min_res), log(max_res), length.out = number_res)))

community_df_res <- furrr::future_map(
  resolutions,
  \(res) {
    set.seed(random_seed)
    community <- igraph::cluster_leiden(
      graph,
      objective_function = 'modularity',
      resolution = res,
      n_iterations = 5L
    )

    modularity <- igraph::modularity(x = graph, membership = community$membership)

    community_df <- data.table::data.table(
      resolution = res,
      node_name = community$names,
      membership = community$membership,
      modularity = modularity
    )
  },
  .progress = .verbose,
  .options = furrr::furrr_options(seed = TRUE)
) %>% data.table::rbindlist(.)

unique(community_df_res[, c("resolution", "modularity")])

community_df_res[, best_modularity := modularity == max(modularity)]


community_df_res[(best_modularity)][, .N, membership] %$% summary(N)
