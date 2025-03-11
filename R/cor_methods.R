# helpers ----

#' Assess the quality of clusters
#'
#' @description Helper function that will assess the quality (defined as the
#' median absolute correlation within the cluster) of a given set of clusters.
#'
#' @param community_df A data.table. Contains node_name (name of the node in the
#' graph) and membership (membership of the node to the respective community).
#' @param correlation_res A data.table. Contains feature_a (first node name),
#' feature_b (second node name) and the r_abs (absolute correlation coefficient)
#' columns.
#'
#' @returns A data.table with the followoing columns:
#'
#' @export
cor_cluster_quality <- function(community_df, correlation_res) {
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

  # Ensure that the keys exist
  data.table::setkey(correlation_res, feature_a, feature_b)

  # Pre-allocate
  cluster_quality <- vector(mode = "list", length = length(community_df_list))

  for (i in seq_along(cluster_genes)) {
    nodes <- cluster_genes[[i]]
    community <- names(cluster_genes)[i]

    node_combinations <- data.table::CJ(feature_a = nodes, feature_b = nodes)

    # Use binary join instead of filtering (much faster for large data)
    # F--k knows why ... ?
    subset_data <- correlation_res_red[node_combinations, on = .(feature_a, feature_b), nomatch = 0]

    median_val <- median(subset_data$r_abs, na.rm = TRUE)
    mad_val <- mad(subset_data$r_abs, na.rm = TRUE)
    size <- length(nodes)

    res <- data.table::data.table(
      module_name = community,
      r_median = median_val,
      r_mad = mad_val,
      r_adjusted = median_val * log10(size),
      size = size
    )

    cluster_quality[[i]] <- res
  }

  cluster_quality <- data.table::rbindlist(cluster_quality)

  cluster_quality
}

# methods - simple correlations ----

#' Prepare correlation-based module detection
#'
#' @description
#' This is the generic function for running a single correlation across your
#' genes.
#'
#' @export
cor_module_processing <- S7::new_generic(
  "cor_module_processing",
  "bulk_coexp"
)

#' @name cor_module_processing
#'
#' @description
#' This function will generate necessarily data for correlation-based detection
#' of gene modules. In this case, this will be based on simple correlation
#' between the genes and the relevant data will be added to the `bulk_coexp`
#' class.
#'
#' @usage ...
#'
#' @return `bulk_coexp` with the needed data for subsequent identification of
#' correlation-based co-expression modules.
#' @param correlation_method String. Option of `c("pearson", "spearman")`.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @export
#'
#' @importFrom magrittr `%>%`
#' @importFrom magrittr `%$%`
#' @importFrom zeallot `%<-%`
#' @import data.table
#'
#' @method cor_module_processing bulk_coexp
S7::method(cor_module_processing, bulk_coexp) <- function(
    bulk_coexp,
    correlation_method = c("pearson", "spearman"),
    .verbose = TRUE) {
  # Checks
  checkmate::assertClass(bulk_coexp, "bixverse::bulk_coexp")
  checkmate::assertChoice(correlation_method, c("pearson", "spearman"))
  checkmate::qassert(.verbose, "B1")

  # Function body
  if (purrr::is_empty(S7::prop(bulk_coexp, "processed_data")[['processed_data']])) {
    warning("No pre-processed data found. Defaulting to the raw data")
    target_mat <- S7::prop(bulk_coexp, "raw_data")
  } else {
    target_mat <- S7::prop(bulk_coexp, "processed_data")[['processed_data']]
  }

  spearman <- if (correlation_method == 'pearson') {
    if (.verbose)
      message("Using Pearson correlations.")
    FALSE
  } else {
    if (.verbose)
      message("Using Spearman correlations.")
    TRUE
  }

  # Calculate the upper triangle of correlation matrix
  cor_diagonal <- rs_cor_upper_triangle(target_mat, spearman = spearman, shift = 1L)

  # Save data to memory friendly R6 class
  cor_data <- upper_triangular_cor_mat$new(
    cor_coef = cor_diagonal,
    features = colnames(target_mat),
    shift = 1L
  )

  correlation_params <- list(
    spearman = spearman,
    type = 'simple'
  )

  S7::prop(bulk_coexp, "processed_data")[["correlation_res"]] <- cor_data
  S7::prop(bulk_coexp, "params")[["correlation_params"]] <- correlation_params

  return(bulk_coexp)
}


# methods - graph-based gene module detection ----

#' Identify correlation-based gene modules via graph methods.
#'
#' @description
#' This is the generic function for identifying gene modules via graph-based
#' methods.
#'
#' @export
cor_module_resolutions <- S7::new_generic(
  "cor_module_resolutions",
  "bulk_coexp"
)



#' @name cor_module_resolutions
#'
#' @param bulk_coexp description
#' @param kernel_bandwidth description
#' @param min_res description
#' @param max_res description
#' @param by description
#' @param random_seed description
#' @param .verbose description
#'
#' @description
#' ...
#'
#' @usage ...
#'
#' @return `bulk_coexp` class. You need to run either
#' [bixverse::cor_module_processing()] or ... .
#' @param correlation_method ...
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @export
#'
#' @importFrom magrittr `%>%`
#' @importFrom magrittr `%$%`
#' @importFrom zeallot `%<-%`
#' @import data.table
#'
#' @method cor_module_resolutions bulk_coexp
S7::method(cor_module_resolutions, bulk_coexp) <- function(bulk_coexp,
                                                           kernel_bandwidth = 0.25,
                                                           min_res = 0.5,
                                                           max_res = 5,
                                                           by = 0.5,
                                                           random_seed = 123L,
                                                           cores = as.integer(parallel::detectCores() / 2),
                                                           parallel = TRUE,
                                                           .verbose = TRUE) {
  # Checks
  checkmate::assertClass(bulk_coexp, "bixverse::bulk_coexp")
  checkmate::qassert(kernel_bandwidth, "N1")
  checkmate::qassert(min_res, "N1")
  checkmate::qassert(max_res, "N1")
  checkmate::qassert(by, "N1")
  checkmate::qassert(random_seed, "I1")
  checkmate::qassert(.verbose, "B1")

  # Early return
  if (purrr::is_empty(S7::prop(bulk_coexp, "processed_data")[['correlation_res']])) {
    warning("No correlation results found. Returning class as is.")
    return(bulk_coexp)
  }

  # Return the data.table from the correlation results
  # correlation_res <- S7::prop(bulk_coexp, "processed_data")$correlation_res$get_data_table() %>%
  #   .[, cor_abs := abs(cor)] %>%
  #   .[, dist := 1 - cor_abs] %>%
  #   .[, dist := data.table::fifelse(dist < 0, 0, dist)] %>%
  #   .[, affinitiy := rs_gaussian_affinity_kernel(x = dist, bandwidth = kernel_bandwidth)] %>%
  #   .[affinitiy >= ]

  dist <- 1 - correlation_res$r_abs
  # Deal with float precision errors because of Rust <> R interface
  dist[dist < 0] <- 0
  affinity <- rs_gaussian_affinity_kernel(x = dist, bandwidth = kernel_bandwidth)

  graph_df <- correlation_res[, c("feature_a", "feature_b")] %>%
    .[, weight := affinity] %>%
    data.table::setnames(.,
                         old = c("feature_a", "feature_b"),
                         new = c("from", "to"))

  graph <- igraph::graph_from_data_frame(graph_df, directed = FALSE)
  graph <- igraph::simplify(graph)

  resolutions <- seq(from = min_res, max_res, by = by)

  if (.verbose)
    message(sprintf("Iterating through %i resolutions", length(resolutions)))
  resolution_results <- purrr::map_dfr(resolutions, \(resolution) {
    set.seed(random_seed)

    community <- igraph::cluster_leiden(graph,
                                        objective_function = 'modularity',
                                        resolution = resolution)

    community_df <- data.table(node_name = community$names,
                               membership = community$membership)

    qc <- cor_cluster_quality(community_df = community_df, correlation_res = correlation_res)

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

  S7::prop(bulk_coexp, "outputs")[['cluser_quality']] <- resolution_results

  return(bulk_coexp)
}

# methods - plotting ----
