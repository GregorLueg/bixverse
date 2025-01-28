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

test_class <- diffuse_seed_genes(test_class, diffusion_vector, 'max')

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

? rs_fast_auc

rs_fast_auc(pos.scores, neg.scores, 10000, 42)

? rs_create_random_aucs

? tied_diffusion

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

# Write the method

network_diffusions = test_class
diffusion_threshold = .5
max_genes = 300L
min_genes = 10L
min_seed_genes = 0L
intial_res = 0.5
seed = 42L
.verbose = T
.max_iters = 100L

# Checks
checkmate::assertClass(network_diffusions, "BIXverse::network_diffusions")
checkmate::qassert(diffusion_threshold, "R1[0,1]")
checkmate::qassert(max_genes, "I1")
checkmate::qassert(min_genes, "I1")
checkmate::qassert(min_seed_genes, "I1")
checkmate::qassert(intial_res, "R1")
checkmate::qassert(seed, "I1")
checkmate::qassert(.verbose, "B1")
checkmate::qassert(.max_iters, "I1")


# Body
## Reduce the graph
diffusion_score <- S7::prop(network_diffusions, "diffusion_res")
if (length(diffusion_score) == 0) {
  warning(
    "The diffusion score has length 0. Likely you did not run the diffusion methods. Returning NULL."
  )
  return(NULL)
}
nodes_to_include <- diffusion_score %>%
  sort(decreasing = T) %>%
  .[1:ceiling(diffusion_threshold * length(diffusion_score))]


red_graph <- igraph::subgraph(S7::prop(network_diffusions, "graph"),
                              names(nodes_to_include))

## First clustering
set.seed(seed)

if (S7::prop(network_diffusions, "params")$directed_graph) {
  red_graph = igraph::as_undirected(red_graph, mode = 'each')
}

first_clusters <- igraph::cluster_leiden(red_graph, resolution = intial_res)

clusters_df <- data.table::data.table(nodeID = first_clusters$names, clusterID = first_clusters$membership)

node_frequency <- clusters_df[, .N, .(clusterID)]

## Subclustering
clusters_with_too_many_nodes <- node_frequency[N > max_genes, clusterID]
final_clusters <- clusters_df[!clusterID %in% clusters_with_too_many_nodes]

for (cluster_i in clusters_with_too_many_nodes) {
  nodes_in_cluster <- clusters_df[clusterID == cluster_i, nodeID]
  finalised_clusters <- data.table()
  # Loop through, until all clusters are below the minimum genes or max iterations is hit
  l <- 1
  while (length(nodes_in_cluster) != 0) {
    set.seed(random_seed + l)
    if (verbose)
      message("Cluster ", i, " gets subclustered. Iter: ", l)
    red_graph_l <- igraph::subgraph(red_graph,
                                    data.table::chmatch(nodes_in_cluster, igraph::V(red_graph)$name))

    clusters_red <- igraph::cluster_leiden(red_graph_l, resolution = intial_res + l * 0.05)

    subclusters <- data.table(nodeID = clusters_red$names,
                              clusterID = clusters_red$membership)
    subclusters_frequency <- subclusters[, .N, .(clusterID)]
    clusters_small_enough <- subclusters_frequency[N <= max_genes, clusterID]

    good_clusters <- subclusters[clusterID %in% clusters_small_enough] %>%
      dplyr::mutate(clusterID = paste0(i, paste(rep("sub", l), collapse = ""), clusterID))

    finalised_clusters <- rbind(finalised_clusters, good_clusters)

    l <- l + 1
    if (l == .max_iters)
      break

    nodes_in_cluster <- setdiff(genes_in_cluster, good_clusters$nodeID)
  }

  final_clusters <- rbind(final_clusters, finalised_clusters)
}

## Add the seed node information based on diffusion type
diffusion_type <- S7::prop(network_diffusions, "params")$diffusion_type

final_node_frequency <- if (diffusion_type == 'single') {
  seed_genes <- S7::prop(network_diffusions, "params")$seed_genes

  final_clusters[, .(cluster_size = length(nodeID),
                     seed_nodes = sum(nodeID %in% seed_genes)), .(clusterID)]
} else {
  seed_genes_set_1 <- S7::prop(network_diffusions, "params")$seed_genes$set_1
  seed_genes_set_2 <- S7::prop(network_diffusions, "params")$seed_genes$set_2

  final_clusters[, .(
    cluster_size = length(nodeID),
    seed_nodes_1 = sum(nodeID %in% seed_genes_set_1),
    seed_nodes_2 = sum(nodeID %in% seed_genes_set_2)
  ), .(clusterID)]
}

## Finalise the clusters
clusters_to_take <- if (diffusion_type == 'single') {
  final_node_frequency[cluster_size >= min_genes &
                         seed_nodes >= min_seed_genes, clusterID]
} else {
  final_node_frequency[cluster_size >= min_genes &
                         seed_nodes_1 >= min_seed_genes &
                         seed_nodes_2 >= min_seed_genes, clusterID]
}

finalised_clusters.clean <- final_clusters[clusterID %in% clusters_to_take]

ks_vals = vector(mode = "numeric", length = length(clusters_to_take))

for (i in seq_along(clusters_to_take)) {
  cluster <- clusters_to_take[i]
  cluster_nodes <- finalised_clusters.clean[clusterID == cluster, nodeID]
  ks <- ks.test(diffusion_score[cluster_nodes],
                diffusion_score[which(!names(diffusion_score) %in% cluster_nodes)],
                alternative = "less")
  ks_vals[i] <- ks$p.value
}

ks_val_df <- data.table(
  clusterID = clusters_to_take,
  ks_pval = ks_vals
)

final_result <- purrr::reduce(list(finalised_clusters.clean, ks_val_df, final_node_frequency),
                              merge,
                              by = 'clusterID') %>%
  .[, `:=`(diffusion_score = diffusion_score[nodeID],
           fdr = p.adjust(ks_pval, method = "fdr"))]

