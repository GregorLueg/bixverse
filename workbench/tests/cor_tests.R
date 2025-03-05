library(devtools)
library(ggplot2)
library(magrittr)

devtools::document()
devtools::load_all()
rextendr::document()

devtools::install()

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

library(recount3)
library(limma)
library(edgeR)

gtex_brain <- recount3::create_rse_manual(
  project = "BRAIN",
  project_home = "data_sources/gtex",
  organism = "human",
  annotation = "gencode_v29",
  type = "gene"
)

coldata <- colData(gtex_brain) |> as.data.frame()
rowdata <- rowData(gtex_brain) |> as.data.frame()

d <- DGEList(assay(gtex_brain))
d <- calcNormFactors(d, method = 'upperquartile')
to_keep <- suppressWarnings(filterByExpr(d))
d <- d[to_keep, ]
d <- cpm(d, log = TRUE)

d <- as.matrix(d)

new_meta_data <- data.table(
  sample_id = rownames(coldata),
  case_control = 'case'
)
d[1:5, 1:5]

cor_test = bulk_coexp(t(d), new_meta_data)

cor_test = preprocess_bulk_coexp(cor_test, hvg = 0.5)

tictoc::tic()
cor_test = cor_module_processing(cor_test, correlation_method = 'spearman')
tictoc::toc()

cor_test = cor_module_identification(cor_test)

kernel_bandwidth = 0.25
min_res = 0.5
max_res = 5
by = 0.5
random_seed = 123L
min_affinity = 0.1
min_genes = 10L
.verbose = TRUE

correlation_res <- S7::prop(cor_test, "processed_data")[['correlation_res']]

dist <- 1 - correlation_res$r_abs
# Deal with float precision errors because of Rust <> R interface
dist[dist < 0] <- 0
affinity <- rs_gaussian_affinity_kernel(x = dist, bandwidth = kernel_bandwidth)

to_keep = (affinity > min_affinity)

graph_df <- correlation_res[, c("feature_a", "feature_b")] %>%
  .[, weight := affinity] %>%
  data.table::setnames(.,
                       old = c("feature_a", "feature_b"),
                       new = c("from", "to")) %>%
  .[(to_keep)]


graph <- igraph::graph_from_data_frame(graph_df, directed = FALSE)
graph <- igraph::simplify(graph)

resolutions <- seq(from = min_res, max_res, by = by)

if (.verbose) message(sprintf("Iterating through %i resolutions", length(resolutions)))

correlation_res_red <- correlation_res[(to_keep)]
data.table::setkey(correlation_res_red, feature_a, feature_b)

future::plan(future::multisession(workers = parallel::detectCores() / 2))

community_df_res <- furrr::future_map(
  resolutions,
  \(res) {
    set.seed(123)
    community <- igraph::cluster_leiden(graph,
                                        objective_function = 'modularity',
                                        resolution = res)

    community_df <- data.table::data.table(
      resolution = res,
      node_name = community$names,
      membership = community$membership
    )
  },
  .progress = TRUE,
  .options = furrr::furrr_options(seed = TRUE)
) %>% data.table::rbindlist(.)

future::plan(future::sequential()); gc()

community_df_res[, combined_id := sprintf("id_%s_%s", resolution, membership)]

clusters_to_keep <- community_df_res[, .N, combined_id][N >= min_genes, combined_id]

community_df_res <- community_df_res[combined_id %in% clusters_to_keep]

cluster_genes <- split(community_df_res$node_name, community_df_res$combined_id)

tictoc::tic()
cluster_quality <- purrr::imap_dfr(cluster_genes, \(nodes, community) {
  c(median, mad) %<-% correlation_res_red[feature_a %in% nodes &
                                            feature_b %in% nodes, c(median(r_abs, na.rm = TRUE), mad(r_abs, na.rm = TRUE))]
  size <- length(nodes)

  res <- data.table::data.table(
    module_name = community,
    r_median = median,
    r_mad = mad,
    r_adjusted = median * log1p(size),
    size = size
  )

  res
})
tictoc::toc()

?CJ



tictoc::tic()
for (i in seq_along(cluster_genes)) {
  nodes <- cluster_genes[[i]]
  community <- names(cluster_genes)[i]

  # Create a temporary data.table of node combinations within this cluster
  node_combinations <- CJ(feature_a = nodes, feature_b = nodes)

  # Use binary join instead of filtering (much faster for large data)
  subset_data <- correlation_res_red[node_combinations, on = .(feature_a, feature_b), nomatch = 0]

  # Calculate statistics
  median_val <- median(subset_data$r_abs, na.rm = TRUE)
  mad_val <- mad(subset_data$r_abs, na.rm = TRUE)
  size <- length(nodes)

  # Add to results
  cluster_quality <- rbindlist(list(
    cluster_quality,
    data.table::data.table(
      module_name = community,
      r_median = median_val,
      r_mad = mad_val,
      r_adjusted = median_val * log1p(size),
      size = size
    )
  ))
}
tictoc::toc()

future::plan(future::multisession(workers = parallel::detectCores() / 2))

tictoc::tic()
# Process clusters in parallel
results_list <- furrr::future_map(seq_along(cluster_genes), function(i) {
  nodes <- cluster_genes[[i]]
  community <- names(cluster_genes)[i]

  # Create indices for filtering
  node_idx <- which(correlation_res_red$feature_a %in% nodes &
                      correlation_res_red$feature_b %in% nodes)

  # Subset using indices (faster than logical filtering)
  r_abs_values <- correlation_res_red$r_abs[node_idx]

  # Calculate statistics
  median_val <- median(r_abs_values, na.rm = TRUE)
  mad_val <- mad(r_abs_values, na.rm = TRUE)
  size <- length(nodes)

  # Return as a single row data.table
  data.table(
    module_name = community,
    r_median = median_val,
    r_mad = mad_val,
    r_adjusted = median_val * log1p(size),
    size = size
  )
})
tictoc::toc()

# Combine results efficiently
cluster_quality <- rbindlist(results_list)
