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

new_meta_data <- data.table::data.table(
  sample_id = rownames(coldata),
  case_control = 'case',
  gtex_subgrp = coldata$gtex.smtsd
)

samples_to_keep <- new_meta_data[gtex_subgrp == "Brain - Hippocampus", sample_id]

cor_test = bulk_coexp(t(d)[samples_to_keep, ], new_meta_data[gtex_subgrp == "Brain - Hippocampus"])

cor_test = preprocess_bulk_coexp(cor_test, mad_threshold = 1)

cor_test = cor_module_processing(cor_test, correlation_method = 'spearman')

cor_test = cor_module_check_res(cor_test)

plot_df <- cor_test@processed_data$feature_meta

ggplot(data = plot_df,
       mapping =  aes(x = MAD)) +
  geom_histogram(mapping = aes(fill = hvg), bins = 100) +
  theme_minimal()

?geom_histogram

kernel_bandwidth = 0.2
min_res = 0.1
max_res = 10
number_res = 15
random_seed = 123L
min_affinity = 0.001
min_genes = 10L
parallel = TRUE
max_workers = as.integer(parallel::detectCores() / 2)
.verbose = TRUE

x <- cor_test

S7::prop(x, "params")[["detection_method"]] == "correlation-based"


graph <- get_cor_graph_single(
  cor_test,
  kernel_bandwidth = kernel_bandwidth,
  min_affinity = min_affinity,
  .verbose = TRUE
)


resolutions <- exp(seq(
  log(min_res),
  log(max_res),
  length.out = number_res
))

if (.verbose) message(sprintf("Iterating through %i resolutions", length(resolutions)))

if(parallel) {
  if (.verbose)
    message(sprintf("Using parallel computation over %i cores.", max_workers))
  future::plan(future::multisession(workers = max_workers))
} else {
  if (.verbose)
    message("Using sequential computation.")
  future::plan(future::sequential())
}

community_df_res <- furrr::future_map(
  resolutions,
  \(res) {
    set.seed(123)
    community <- igraph::cluster_leiden(graph,
                                        objective_function = 'modularity',
                                        resolution = res,
                                        n_iterations = 5L)

    modularity <- igraph::modularity(x = graph, membership = community$membership)

    community_df <- data.table::data.table(
      resolution = res,
      node_name = community$names,
      membership = community$membership,
      modularity = modularity
    )
  },
  .progress = TRUE,
  .options = furrr::furrr_options(seed = TRUE)
) %>% data.table::rbindlist(.)

future::plan(future::sequential()); gc()

community_df_res[, combined_id := sprintf("id_%s_%s", resolution, membership)]

cluster_summary <- community_df_res[, .N, combined_id] %>%
  .[, good_clusters := N >= min_genes] %>%
  data.table::merge.data.table(.,
                               unique(community_df_res[, c('resolution', 'combined_id')]),
                               by.x = 'combined_id',
                               by.y = 'combined_id',
  ) %>%
  .[, .(
    good_clusters = sum(good_clusters),
    avg_size = mean(N),
    max_size = max(N)
  ), resolution]

resolution_results <- community_df_res[, .(no_clusters = length(unique(membership)),
                                           modularity = unique(modularity)), resolution] %>%
  data.table::merge.data.table(., cluster_summary, by.x = 'resolution', by.y = 'resolution')





clusters_to_keep <- community_df_res[, .N, combined_id][N >= min_genes, combined_id]

x <- community_df_res[, .N, .(resolution, membership)]

community_df_res <- community_df_res[combined_id %in% clusters_to_keep]

community_df_res[, c("resolution", "modularity")] %>% unique()

community_df_res[, .(no_clusters = length(unique(membership)),
                     modularity = unique(modularity)), resolution]

length(unique(community_df_res$node_name))

# Into function

tictoc::tic()
cluster_genes <- split(community_df_res$node_name, community_df_res$combined_id)

data.table::setkey(correlation_res, feature_a, feature_b)

cluster_quality <- vector(mode = "list", length = length(cluster_genes))

for (i in seq_along(cluster_genes)) {
  nodes <- cluster_genes[[i]]
  community <- names(cluster_genes)[i]
  node_combinations <- data.table::CJ(feature_a = nodes, feature_b = nodes)
  subset_data <- correlation_res[node_combinations, on = .(feature_a, feature_b), nomatch = 0]

  # subset_data <- correlation_res[feature_a %in% nodes & feature_b %in% nodes]

  mean_abs_cor <- median(subset_data$cor_abs)
  sd_abs_cor <- mad(subset_data$cor_abs)
  size <- length(nodes)

  results <- data.table::data.table(
    module_name = community,
    mean_abs_cor = mean_abs_cor,
    sd_abs_cor = sd_abs_cor,
    size = size
  )

  cluster_quality[[i]] <- results
}

cluster_quality <- data.table::rbindlist(cluster_quality)
tictoc::toc()

cluster_quality_final = merge(
  unique(community_df_res[, c("resolution", "membership", "combined_id", "modularity")]),
  cluster_quality,
  by.x = 'combined_id', by.y = 'module_name')

head(cluster_quality_final)

cluster_quality_summarised <- cluster_quality_final[, .(
  mom_abs_cor = median(mean_abs_cor),
  mean_of_sd_abs_cor = median(sd_abs_cor),
  mean_size = mean(size),
  modularity = unique(modularity),
  cluster_no = .N
), resolution] %>% data.table::setorder(., resolution)

head(cluster_quality_summarised)

ggplot(
  data = cluster_quality_summarised,
  mapping = aes(
    x = resolution,
    y = modularity
  )
) +
  geom_point(aes(color = mom_abs_cor, size = sqrt(cluster_no)))



i <- 4

nodes <- cluster_genes[[i]]
community <- names(cluster_genes)[i]

tictoc::tic()
node_combinations <- data.table::CJ(feature_a = nodes, feature_b = nodes)
subset_data <- correlation_res[node_combinations, on = .(feature_a, feature_b), nomatch = 0]
tictoc::toc()

tictoc::tic()

tictoc::toc()


quantile(subset_data$cor_abs, 0.75)[1]

