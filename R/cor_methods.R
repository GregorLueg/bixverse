# methods - helpers ------------------------------------------------------------

## graph generation ------------------------------------------------------------

#' @title Get correlation-based graph
#'
#' @description
#' Helper function to get a correlation-based igraph from the class
#'
#' @param `bulk_coexp` The class, see [bixverse::bulk_coexp()].
#' @param kernel_bandwidth Numerical. The bandwidth to use for the affinity
#' kernel
#' @param min_affinity Numerical. Minimum affinity needed to keep the edge.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @return igraph given the provided parameters.
get_cor_graph_single <- S7::new_generic(
  name = 'get_cor_graph_single',
  dispatch_args = 'bulk_coexp',
  fun = function(bulk_coexp,
                 kernel_bandwidth,
                 min_affinity,
                 .verbose) {
    S7::S7_dispatch()
  }
)

#' @method get_cor_graph_single bulk_coexp
S7::method(get_cor_graph_single, bulk_coexp) <- function(bulk_coexp,
                                                         kernel_bandwidth,
                                                         min_affinity,
                                                         .verbose) {
  # Checks
  checkmate::assertClass(bulk_coexp, "bixverse::bulk_coexp")
  checkmate::qassert(kernel_bandwidth, "R1[0,1]")
  checkmate::qassert(min_affinity, "R1[0,1]")
  checkmate::qassert(.verbose, "B1")
  # Function body
  cor_res <- S7::prop(bulk_coexp, "processed_data")$correlation_res
  graph_df <- cor_res$get_data_table(.verbose = .verbose) %>%
    .[, cor_abs := abs(cor)] %>%
    .[, dist := 1 - cor_abs] %>%
    .[, dist := data.table::fifelse(dist < 0, 0, dist)] %>%
    .[, affinity := rs_gaussian_affinity_kernel(x = dist, bandwidth = kernel_bandwidth)] %>%
    .[affinity >= min_affinity] %>%
    .[, c("feature_a", "feature_b", "affinity")] %>%
    data.table::setnames(
      .,
      old = c("feature_a", "feature_b", "affinity"),
      new = c("from", "to", "weight")
    )

  if (.verbose)
    message("Generating correlation-based graph.")

  graph <- igraph::graph_from_data_frame(graph_df, directed = FALSE)
  graph <- igraph::simplify(graph)

  graph
}

## getters ---------------------------------------------------------------------

#' @title Return the resolution results
#'
#' @description
#' Getter function to get the resolution results (if available).
#'
#' @param `bulk_coexp` The class, see [bixverse::bulk_coexp()].
#'
#' @return If resolution results were found, returns the data.table. Otherwise,
#' throws a warning and returns NULL.
#'
#' @export
get_resolution_res <- S7::new_generic(
  name = 'get_resolution_res',
  dispatch_args = 'bulk_coexp',
  fun = function(bulk_coexp) {
    S7::S7_dispatch()
  }
)

#' @export
#' @method get_resolution_res bulk_coexp
S7::method(get_resolution_res, bulk_coexp) <- function(bulk_coexp) {
  # Checks
  checkmate::assertClass(bulk_coexp, "bixverse::bulk_coexp")
  # Body
  resolution_results <- S7::prop(bulk_coexp, "outputs")[['resolution_results']]
  if (is.null(resolution_results))
    warning("No resolution results found. Did you run cor_module_check_res()? Returning NULL.")

  resolution_results
}

## plotting --------------------------------------------------------------------

#' @title Plot the resolution results.
#'
#' @description
#' Plots the resolution results (if they can be found in the class). The x-axis
#' reflects the
#'
#' @param `bulk_coexp` The class, see [bixverse::bulk_coexp()].
#'
#' @return If resolution results were found, returns the data.table. Otherwise,
#' throws a warning and returns NULL.
#'
#' @export
plot_resolution_res <- S7::new_generic(
  name = 'plot_resolution_res',
  dispatch_args = 'bulk_coexp',
  fun = function(bulk_coexp, print_head = TRUE) {
    S7::S7_dispatch()
  }
)

#' @export
#'
#' @import ggplot2
#'
#' @method plot_resolution_res bulk_coexp
S7::method(plot_resolution_res, bulk_coexp) <- function(bulk_coexp, print_head = TRUE) {
  # Checks
  checkmate::assertClass(bulk_coexp, "bixverse::bulk_coexp")
  checkmate::qassert(print_head, "B1")
  # Body
  plot_df <- S7::prop(bulk_coexp, "outputs")[['resolution_results']]
  if (is.null(plot_df)) {
    warning("No resolution results found. Did you run cor_module_check_res()? Returning NULL.")
    return(NULL)
  }

  plot_df <- data.table::setorder(plot_df, -modularity)

  # if (print_head)
  #   print(head(plot_df))

  p <- ggplot(data = plot_df,
              mapping =  aes(x = resolution, y = modularity)) +
    geom_point(
      mapping = aes(size = log10(good_clusters), fill = log10(avg_size)),
      shape = 21,
      alpha = .7
    ) +
    xlab("Leiden cluster resolution") +
    ylab("Modularity") +
    theme_minimal() +
    scale_fill_viridis_c() +
    scale_size_continuous(range = c(2, 8)) +
    labs(size = "Number of good clusters (log10)", fill = "Average cluster size (log10)") +
    ggtitle("Resolution vs. modularity", subtitle = 'With cluster number and size')

  p
}

# methods - simple correlations ------------------------------------------------

#' @title Prepare correlation-based module detection
#'
#' @description
#' This function will calculate the correlation coefficients between the genes,
#' using the highly variable genes (if available, otherwise the function will
#' use the raw data). The data will be stored in a memory-efficient format
#' in the properties of the class.
#'
#' @return `bulk_coexp` The class, see [bixverse::bulk_coexp()]. Ideally, you
#' should run [bixverse::preprocess_bulk_coexp()] before applying this class.
#' @param correlation_method String. Option of `c("pearson", "spearman")`.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @return The class with added data to the properties for subsequent functions.
#'
#' @export
cor_module_processing <- S7::new_generic(
  name = "cor_module_processing",
  dispatch_args = "bulk_coexp",
  fun = function(bulk_coexp,
                 correlation_method = c("pearson", "spearman"),
                 .verbose = TRUE) {
    S7::S7_dispatch()
  }
)


#' @export
#'
#' @importFrom magrittr `%>%`
#' @importFrom magrittr `%$%`
#' @importFrom zeallot `%<-%`
#' @import data.table
#'
#' @method cor_module_processing bulk_coexp
S7::method(cor_module_processing, bulk_coexp) <- function(bulk_coexp,
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

  correlation_params <- list(spearman = spearman, type = 'simple')

  S7::prop(bulk_coexp, "processed_data")[["correlation_res"]] <- cor_data
  S7::prop(bulk_coexp, "params")[["correlation_params"]] <- correlation_params
  S7::prop(bulk_coexp, "params")["detection_method"] <- "correlation-based"

  return(bulk_coexp)
}


# methods - graph-based gene module detection ----------------------------------

#' @title Iterate through Leiden resolutions for graph-based community detection.
#'
#' @description
#' This function will identify gene modules based on affinity graphs from the
#' single correlation or differential correlation methods. Briefly, in the case
#' of single correlation, the graph is generated based on the absolute
#' correlation coefficients that are subjected to a Gaussian affinity kernel.
#' TODO: Write what is being done for differential correlation methods.
#' This reduces spurious correlations and leaves a sparsely connected graph.
#' Subsequently, Leiden community detection is applied through a range of
#' resolutions that the user can define. The function then returns meta
#' information about the resolutions (which can also be plotted) to identify
#' the best suitable resolution parameter to identify co-expression modules.
#'
#' @param `bulk_coexp` The class, see [bixverse::bulk_coexp()].
#' @param min_res Numeric. The minimum resolution to test for the Leiden
#' community detection.
#' @param max_res Numeric. The maximum resolution to test for the Leiden
#' community detection.
#' @param number_res Integer. Number of resolutions to test.
#' The resolutions will be spread in a logarithmic fashion over `min_res` to
#' `max_res`.
#' @param random_seed Integer. Random seed.
#' @param kernel_bandwidth Numeric. !This parameter is only relevant for simple
#' correlation modules!. The bandwidth for the affinity kernel. Needs to be
#' value between 0 and 1.
#' @param min_affinity Float. !This parameter is only relevant for simple
#' correlation modules! This parameter will remove edges below this affinity
#' threshold. Needs to be a value between 0 and 1.
#' @param min_genes Integer. Minimum number of genes that should be in a
#' community.
#' @param parallel Boolean. Parallelise the Leiden clustering.
#' @param max_workers Integer. Maximum number of workers to use if parallel is
#' set to `TRUE`.
#' @param .verbose Controls the verbosity of the function.
#'
#' @return The class with added data to the properties.
#'
#' @export
cor_module_check_res <- S7::new_generic(
  name = "cor_module_check_res",
  dispatch_args = "bulk_coexp",
  fun = function(bulk_coexp,
                 min_res = 0.1,
                 max_res = 10,
                 number_res = 15L,
                 random_seed = 123L,
                 kernel_bandwidth = 0.2,
                 min_affinity = 0.001,
                 min_genes = 10L,
                 parallel = TRUE,
                 max_workers = as.integer(parallel::detectCores() / 2),
                 .verbose = TRUE) {
    S7::S7_dispatch()
  }
)


#' @export
#'
#' @importFrom magrittr `%>%`
#' @import data.table
#'
#' @method cor_module_check_res bulk_coexp
S7::method(cor_module_check_res, bulk_coexp) <- function(bulk_coexp,
                                                         min_res = 0.1,
                                                         max_res = 10,
                                                         number_res = 15L,
                                                         random_seed = 123L,
                                                         kernel_bandwidth = 0.2,
                                                         min_affinity = 0.001,
                                                         min_genes = 10L,
                                                         parallel = TRUE,
                                                         max_workers = as.integer(parallel::detectCores() / 2),
                                                         .verbose = TRUE) {
  # Checks
  checkmate::assertClass(bulk_coexp, "bixverse::bulk_coexp")
  checkmate::qassert(min_res, "N1")
  checkmate::qassert(max_res, "N1")
  checkmate::qassert(random_seed, "I1")
  checkmate::qassert(kernel_bandwidth, "R1[0,1]")
  checkmate::qassert(min_affinity, "R1[0,1]")
  checkmate::qassert(number_res, "I1")
  checkmate::qassert(min_genes, "I1")
  checkmate::qassert(parallel, "B1")
  checkmate::qassert(max_workers, "I1")
  checkmate::qassert(.verbose, "B1")

  # Early return
  if (purrr::is_empty(S7::prop(bulk_coexp, "processed_data")[['correlation_res']])) {
    warning("No correlation results found. Returning class as is.")
    return(bulk_coexp)
  }

  graph <- if (S7::prop(bulk_coexp, "params")[["detection_method"]] == "correlation-based") {
    get_cor_graph_single(
      bulk_coexp,
      kernel_bandwidth = kernel_bandwidth,
      min_affinity = min_affinity,
      .verbose = .verbose
    )
  }

  resolutions <- exp(seq(log(min_res), log(max_res), length.out = number_res))

  if (.verbose)
    message(sprintf("Iterating through %i resolutions", length(resolutions)))

  if (parallel) {
    if (.verbose)
      message(sprintf("Using parallel computation over %i cores.", max_workers))

    # future plan funkiness
    assign(".temp_workers", max_workers, envir = .GlobalEnv)
    on.exit(rm(".temp_workers", envir = .GlobalEnv))

    future::plan(future::multisession(workers = .temp_workers))
  } else {
    if (.verbose)
      message("Using sequential computation.")
    future::plan(future::sequential())
  }

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
  # To make the message trace prettier
  if(.verbose)
    message("")

  future::plan(future::sequential())
  gc()

  community_df_res[, combined_id := sprintf("id_%s_%s", resolution, membership)]

  cluster_summary <- community_df_res[, .N, combined_id] %>%
    .[, good_clusters := N >= min_genes] %>%
    data.table::merge.data.table(.,
                                 unique(community_df_res[, c('resolution', 'combined_id')]),
                                 by.x = 'combined_id',
                                 by.y = 'combined_id',) %>%
    .[, .(
      good_clusters = sum(good_clusters),
      avg_size = mean(N),
      max_size = max(N)
    ), resolution]

  resolution_results <- community_df_res[, .(no_clusters = length(unique(membership)),
                                             modularity = unique(modularity)), resolution] %>%
    data.table::merge.data.table(., cluster_summary, by.x = 'resolution', by.y = 'resolution')

  # Store the parameters
  graph_params <- list(
    kernel_bandwidth = kernel_bandwidth,
    min_affinity = min_affinity,
    no_nodes = length(igraph::V(graph)),
    no_edges = length(igraph::E(graph))
  )

  # Assign stuff
  S7::prop(bulk_coexp, "params")[["correlation_graph"]] <- graph_params
  S7::prop(bulk_coexp, "outputs")[['resolution_results']] <- resolution_results
  S7::prop(bulk_coexp, "outputs")[['cor_graph']] <- graph

  return(bulk_coexp)
}

#' @title Identify correlation-based gene modules via graphs.
#'
#' @description
#' This function leverages graph-based clustering to identify gene co-expression
#' modules. The class has the option to sub-cluster large communities within
#' their respective sub graphs, akin to the approach taken by Barrio-Hernandez,
#' et al.
#'
#' @param `bulk_coexp` The class, see [bixverse::bulk_coexp()].
#' @param resolution The Leiden resolution parameter you wish to use. If NULL,
#' it will use the optimal one identified by [bixverse::cor_module_check_res()].
#' If nothing can be found, will default to 1.
#' @param min_size Integer. Minimum size of the communities.
#' @param max_size Integer. Maximum size of the communities.
#' @param subclustering Boolean. Shall after a first clustering communities that
#' are too large be further sub clustered. Defaults to `TRUE`.
#' @param random_seed Integer. Random seed.
#' @param .kernel_bandwidth Numeric. A value between 0 and 1. !This parameter
#' only plays a role if you have not run [bixverse::cor_module_check_res()]!.
#' @param .min_affinity Numeric. A value between 0 and 1. !This parameter only
#' plays a role if you have not run [bixverse::cor_module_check_res()]!.
#' @param .max_iters Integer. If sub clustering is set to `TRUE`, what shall be the
#' maximum number of iterations. Defaults to 100L.
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @return The class with added data to the properties.
#'
#' @references Barrio-Hernandez, et al., Nat Genet, 2023.
#'
#' @export
cor_module_final_modules <- S7::new_generic(
  name = "cor_module_final_modules",
  dispatch_args = "bulk_coexp",
  fun = function(bulk_coexp,
                 resolution = NULL,
                 min_size = 10L,
                 max_size = 500L,
                 subclustering = TRUE,
                 random_seed = 123L,
                 .kernel_bandwidth = 0.2,
                 .min_affinity = 0.001,
                 .max_iters = 100L,
                 .verbose = TRUE) {
    S7::S7_dispatch()
  }
)

#' @export
#'
#' @importFrom magrittr `%>%`
#' @import data.table
#'
#' @method cor_module_final_modules bulk_coexp
S7::method(cor_module_final_modules, bulk_coexp) <- function(bulk_coexp,
                                                             resolution = NULL,
                                                             min_size = 10L,
                                                             max_size = 500L,
                                                             subclustering = TRUE,
                                                             random_seed = 123L,
                                                             .kernel_bandwidth = 0.2,
                                                             .min_affinity = 0.001,
                                                             .max_iters = 100L,
                                                             .verbose = TRUE) {
  # Checks
  checkmate::assertClass(bulk_coexp, "bixverse::bulk_coexp")
  checkmate::qassert(resolution, c("0", "N1"))
  checkmate::qassert(min_size, "I1")
  checkmate::qassert(max_size, "I1")
  checkmate::qassert(subclustering, "B1")
  checkmate::qassert(random_seed, "I1")
  checkmate::qassert(.kernel_bandwidth, "R1[0,1]")
  checkmate::qassert(.min_affinity, "R1[0,1]")
  checkmate::qassert(.max_iters, "I1")
  checkmate::qassert(.verbose, "B1")

  # Early return
  if (purrr::is_empty(S7::prop(bulk_coexp, "processed_data")[['correlation_res']])) {
    warning("No correlation results found. Returning class as is.")
    return(bulk_coexp)
  }

  # Get the graph
  if (is.null(S7::prop(bulk_coexp, "outputs")[["cor_graph"]])) {
    # Deal with the case a graph was not yet generated...
    warning(
      paste(
        "No correlation graph found. Did you run cor_module_check_res()?",
        "Generating correlation graph based on standard parameters.",
        sep = "\n"
      )
    )
    graph <- if (S7::prop(bulk_coexp, "params")[["detection_method"]] == "correlation-based") {
      get_cor_graph_single(
        bulk_coexp,
        kernel_bandwidth = .kernel_bandwidth,
        min_affinity = .min_affinity,
        .verbose = .verbose
      )
    }

    graph_params <- list(
      kernel_bandwidth = kernel_bandwidth,
      min_affinity = min_affinity,
      no_nodes = length(igraph::V(graph)),
      no_edges = length(igraph::E(graph))
    )

    S7::prop(bulk_coexp, "params")[["correlation_graph"]] <- graph_params

  } else {
    graph <- S7::prop(bulk_coexp, "outputs")[["cor_graph"]]
  }

  # Final resolution
  if (is.null(resolution)) {
    resolution_results <- S7::prop(x, "outputs")[["resolution_results"]]
    final_resolution <- if (!is.null(resolution_results)) {
      if (.verbose)
        message("Using resolution with best modularity.")
      resolution_results[modularity == max(modularity), resolution]
    } else {
      warning("No resolution results found and none provided. Will default to a resolution of 1.")
      1
    }
  }

  # Do a first clustering
  set.seed(random_seed)

  final_gene_communities <- igraph::cluster_leiden(
    graph,
    objective_function = 'modularity',
    resolution = final_resolution,
    n_iterations = 5L
  )

  clusters_df <- data.table::data.table(node_id = final_gene_communities$names,
                                        cluster_id = final_gene_communities$membership)

  node_frequency <- clusters_df[, .N, .(cluster_id)]

  # If sub clustering is active, do that
  if (subclustering) {
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

  # Finalise the data
  final_communities <- node_frequency_updated[N <= max_size &
                                                N >= min_size, cluster_id]
  final_clusters_filtered <- final_clusters[cluster_id %in% final_communities]

  cluster_name_prettifier <- setNames(
    paste("cluster", seq_along(
      unique(final_clusters_filtered$cluster_id)
    ), sep = "_"),
    unique(final_clusters_filtered$cluster_id)
  )

  final_clusters_filtered[, cluster_id := cluster_name_prettifier[cluster_id]]

  S7::prop(bulk_coexp, "final_results") <- final_clusters_filtered

  return(bulk_coexp)
}
