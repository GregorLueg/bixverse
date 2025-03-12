library(devtools)
library(ggplot2)
library(magrittr)

devtools::document()
devtools::load_all()
rextendr::document()
devtools::check()

?S7::new_generic

syn_data = synthetic_signal_matrix()

X = t(syn_data$mat)

meta_data = data.table::data.table(
  sample_id = names(syn_data$group),
  case_control = 'case',
  grp = syn_data$group
)

cor_test = bulk_coexp(X, meta_data)

cor_test = preprocess_bulk_coexp(cor_test)

cor_test = cor_module_processing(cor_test, correlation_method = 'spearman')



cor_test@processed_data$correlation_res$get_data_table()

cor_test = cor_module_identification(cor_test)

cor_test@outputs$cluser_quality

plot_df <- cor_test@outputs$cluser_quality

head(plot_df)

ggplot(data = plot_df,
       mapping = aes(x = res, y = r_median_of_adjust)) +
  geom_point(mapping = aes(size = log10(median_size), fill = r_weighted_median),
             shape = 21) +
  theme_minimal() +
  xlab("Leiden resolution") +
  ylab("Median adjusted r")

plot_df

# Test on real data ----

gtex_brain <- recount3::create_rse_manual(
  project = "BRAIN",
  project_home = "data_sources/gtex",
  organism = "human",
  annotation = "gencode_v29",
  type = "gene"
)

coldata <- SummarizedExperiment::colData(gtex_brain) |> as.data.frame()
rowdata <- SummarizedExperiment::rowData(gtex_brain) |> as.data.frame()

d <- edgeR::DGEList(SummarizedExperiment::assay(gtex_brain))
d <- edgeR::calcNormFactors(d, method = 'upperquartile')
to_keep <- suppressWarnings(edgeR::filterByExpr(d))
d <- d[to_keep, ]
d <- edgeR::cpm(d, log = TRUE)

d <- as.matrix(d)

dim(d)

new_meta_data <- data.table::data.table(
  sample_id = rownames(coldata),
  case_control = 'case',
  gtex_subgrp = coldata$gtex.smtsd
)

samples_to_keep <- new_meta_data[gtex_subgrp == "Brain - Hippocampus", sample_id]

dim(t(d)[samples_to_keep, ])

tictoc::tic()
cor_test = bulk_coexp(t(d)[samples_to_keep, ], new_meta_data[gtex_subgrp == "Brain - Hippocampus"]) %>%
  preprocess_bulk_coexp(., hvg = .5) %>%
  cor_module_processing(., correlation_method = 'spearman') %>%
  cor_module_check_res(.) %>%
  cor_module_final_modules(., min_size = 25)

tictoc::toc()

cor_test <- cor_module_final_modules(cor_test, min_size = 25L)

go_data <- get_go_human_data()

go_data <- gene_ontology_data(go_data, min_genes = 5L)

final_results <- cor_test@final_results
final_results[, node_id := gsub("[.].*", "", node_id)]

final_results[, .N, cluster_id]

final_results_list <- split(final_results$node_id, final_results$cluster_id)

tictoc::tic()
go_results <- bixverse::gse_go_elim_method_list(go_data, final_results_list)
tictoc::toc()


head(go_results)

merge()

plot_resolution_res(cor_test)

table(cor_test@final_results$cluster_id)

devtools::document()
devtools::load_all()

cor_test <- cor_module_check_res(cor_test)

?get_params

get_params(cor_test, TRUE, TRUE)

S7::prop(x, "params")[['correlation_graph']][['no_nodes']]

x <- cor_test
resolution = NULL
min_size = 10L
max_size = 500L
subclustering = TRUE
random_seed = 123L
.kernel_bandwidth = 0.2
.min_affinity = 0.001
.max_iters = 100L
.verbose = TRUE

if (is.null(S7::prop(x, "outputs")[["cor_graph"]])) {
  warning(
    paste(
      "No correlation graph found. Did you run cor_module_check_res()?",
      "Generating correlation graph based on standard parameters.",
      sep = "\n"
    )
  )

  # TODO Need to deal with differential correlations here...

  graph <- if(S7::prop(x, "params")[["detection_method"]] == "correlation-based") {
    get_cor_graph_single(
      bulk_coexp,
      kernel_bandwidth = .kernel_bandwidth,
      min_affinity = .min_affinity,
      .verbose = .verbose
    )
  }


} else {
  graph <- S7::prop(x, "outputs")[["cor_graph"]]
}

if (is.null(resolution)) {
  resolution_results <- S7::prop(x, "outputs")[["resolution_results"]]
  final_resolution <- if (!is.null(resolution_results)) {
    if (.verbose)
      message("Using resolution with best modularity.")
    resolution_results[modularity == max(modularity), resolution]
  } else {
    warning("No resolution results found. Will default to a resolution of 1.")
    1
  }
}

set.seed(random_seed)

final_gene_communities <- igraph::cluster_leiden(
  graph,
  objective_function = 'modularity',
  resolution = final_resolution,
  n_iterations = 5L
)

clusters_df <- data.table::data.table(
  node_id = final_gene_communities$names, cluster_id = final_gene_communities$membership
)

node_frequency <- clusters_df[, .N, .(cluster_id)]

if(subclustering) {
  if (.verbose)
    message("Sub-clustering larger communities until they are below max_size.")
  clusters_with_too_many_nodes <- node_frequency[N > max_size, cluster_id]
  final_clusters <- clusters_df[!cluster_id %in% clusters_with_too_many_nodes]

  for (i in seq_along(clusters_with_too_many_nodes)) {
    cluster_i <- clusters_with_too_many_nodes[i]
    nodes_in_cluster <- clusters_df[cluster_id == cluster_i, node_id]
    finalised_clusters <- data.table()
    # Loop through, until all clusters are below the minimum genes or max
    # iterations is hit
    l <- 1
    while (length(nodes_in_cluster) != 0) {
      set.seed(random_seed + l)

      sub_graph_l <- igraph::subgraph(graph,
                                      data.table::chmatch(nodes_in_cluster, igraph::V(graph)$name))

      # Restarting at a very small resolution
      clusters_red <- igraph::cluster_leiden(sub_graph_l, resolution = 0.1 + l * 0.05)

      subclusters <- data.table(node_id = clusters_red$names,
                                cluster_id = clusters_red$membership)

      subclusters_frequency <- subclusters[, .N, .(cluster_id)]
      clusters_small_enough <- subclusters_frequency[N <= max_size, cluster_id]

      good_clusters <- subclusters[cluster_id %in% clusters_small_enough] %>%
        dplyr::mutate(cluster_id = paste0(i, paste(rep("sub", l), collapse = ""), cluster_id))

      finalised_clusters <- rbind(finalised_clusters, good_clusters)

      l <- l + 1
      if (l == .max_iters) {
        break
      }

      nodes_in_cluster <- setdiff(nodes_in_cluster, good_clusters$node_id)
    }

    final_clusters <- rbind(final_clusters, finalised_clusters)
  }

  node_frequency_updated <- final_clusters[, .N, .(cluster_id)]
} else {
  node_frequency_updated <- node_frequency
  final_clusters <- clusters_df
}

final_communities <- node_frequency_updated[N <= max_size & N >= min_size, cluster_id]
final_clusters_filtered <- final_clusters[cluster_id %in% final_communities]

cluster_name_prettifier <- setNames(
  paste("cluster", seq_along(unique(
    final_clusters_filtered$cluster_id
  )), sep = "_"),
  unique(final_clusters_filtered$cluster_id)
)

final_clusters_filtered[, cluster_id := cluster_name_prettifier[cluster_id]]

data <- t(d)[, 1:15000]
features <- colnames(data)[1:15000]

dim(data)

tictoc::tic()
# Rust
flattened_cor <- rs_cor_upper_triangle(data, shift = 1L, spearman = T)

upper_triangle <- upper_triangular_cor_mat$new(cor_coef = flattened_cor,
                                               features = features,
                                               shift = 1L)

matrix <- upper_triangle$get_cor_matrix()
tictoc::toc()

matrix[1:5, 1:5]

tictoc::tic()
matrix_r = cor(data, method = 'spearman')
tictoc::toc()

matrix_r[1:5, 1:5]
