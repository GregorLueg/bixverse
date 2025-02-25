library(devtools)

devtools::document()
devtools::load_all()
rextendr::document()


syn_data = synthetic_signal_matrix()

X = t(syn_data$mat)

meta_data = data.table::data.table(
  sample_id = names(syn_data$group),
  case_control = 'case',
  grp = syn_data$group
)


cor_test = bulk_coexp(X, meta_data)

cor_test = preprocess_bulk_coexp(cor_test)



# Calculate cluster quality
cluster_quality <- function(community_df, correlation_res) {
  # Checks
  checkmate::assertDataTable(community_df)
  checkmate::assertDataTable(correlation_res)
  checkmate::assertNames(names(community_df),
                         must.include = c("node_name", "membership"))
  checkmate::assertNames(names(correlation_res),
                         must.include = c("feature_a", "feature_b", "r_abs"))
  # Function
  community_df_list <- split(community_df$node_name, community_df$membership)
  community_df_list <- purrr::keep(community_df_list, \(x) length(x) > 1)

  cluster_quality <- purrr::imap_dfr(community_df_list, \(nodes, community) {
    c(median, mad) %<-% correlation_res[feature_a %in% nodes &
                                          feature_b %in% nodes, c(median(r_abs, na.rm = TRUE), mad(r_abs, na.rm = TRUE))]
    size <- length(nodes)
    data.table(
      module_name = community,
      r_median = median,
      r_mad = mad,
      size = size,
      r_adjusted = median * log1p(size)
    )
  })
  cluster_quality
}


tictoc::tic()
spearman <- TRUE

target_mat <- S7::prop(cor_test, "processed_data")[['processed_data']]

with_diagonal <- FALSE

shift <- ifelse(with_diagonal, 0L, 1L)

cor_diagonal <- rs_cor_upper_triangle(target_mat, spearman = spearman, shift = shift)

feature_names <- colnames(target_mat)
total_len <- length(feature_names)

feature_a <- purrr::map(1:total_len, \(idx) {
  res <- if (with_diagonal) {
    rep(feature_names[[idx]], total_len - idx + 1)
  } else {
    rep(feature_names[[idx]], total_len - idx)
  }
  res
})
feature_a <- do.call(c, feature_a)

feature_b <- purrr::map(1:total_len, \(idx) {
  res <- if (with_diagonal) {
    feature_names[idx:total_len]
  } else {
    remaining <- idx + 1
    if (remaining <= total_len) {
      feature_names[(idx + 1):total_len]
    } else {
      NA
    }
  }
  res
})
feature_b <- na.omit(do.call(c, feature_b))

correlation_res <- data.table(
  feature_a = feature_a,
  feature_b = feature_b,
  r = cor_diagonal,
  r_abs = abs(cor_diagonal)
)

dist <- 1 - correlation_res$r_abs
dist[dist < 0] <- 0

rs_gauss_kernel <- rs_gaussian_affinity_kernel(dist, 0.25)


graph_df = correlation_res[, c("feature_a", "feature_b")][, weight := rs_gauss_kernel] %>%
  data.table::setnames(., old = c("feature_a", "feature_b"), new = c("from", "to"))

graph <- igraph::graph_from_data_frame(graph_df, directed = FALSE)
graph <- igraph::simplify(graph)



resolutions <- seq(from = 0.5, to = 5, by = .5)

resolution_results <- purrr::map(resolutions, \(resolution) {

  community <- igraph::cluster_leiden(graph, objective_function = 'modularity', resolution = resolution)

  community_df <- data.table(
    node_name = community$names,
    membership = community$membership
  )

  qc <- cluster_quality(community_df = community_df, correlation_res = correlation_res)

  res <- qc[, .(
    res = resolution,
    n = .N,
    "median_size" = median(size, na.rm = TRUE),
    "r_weighted_median" = sum(r_median * size) / sum(size),
    "r_median_of_medians" = median(r_median, na.rm = TRUE),
    "r_median_of_adjust" = median(r_adjusted, na.rm = TRUE)
  )]

  res
})

resolution_results <- rbindlist(resolution_results)

plot(resolution_results$res, resolution_results$r_median_of_adjust,
     xlab = "Leiden resolution",
     ylab = "Median adjusted r")

