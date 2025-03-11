# methods - helpers ------------------------------------------------------------

## graph generation ------------------------------------------------------------

#' Get correlation-based graph
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

  if(.verbose)
    message("Generating correlation-based graph.")

  graph <- igraph::graph_from_data_frame(graph_df, directed = FALSE)
  graph <- igraph::simplify(graph)

  graph
}

## getters ---------------------------------------------------------------------

#' Get the resolution results
#'
#' @description
#' Getter function to get the
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
  if(is.null(resolution_results))
    warning("No resolution results found. Did you run cor_module_check_res()? Returning NULL.")

  resolution_results
}

## plotting --------------------------------------------------------------------

#' Get the resolution results
#'
#' @description
#' Helper function to get a correlation-based igraph from the class
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
  fun = function(bulk_coexp,
                 print_head = TRUE) {
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

#' Prepare correlation-based module detection
#'
#' @description
#' This function will calculate the correlation coefficients between the genes,
#' using the highly variable genes (if available, otherwise the function will
#' use the raw data). The data will be stored in a memory-efficient format as a
#' vector within the class.
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
  S7::prop(bulk_coexp, "params")["detection_method"] <- "correlation-based"

  return(bulk_coexp)
}


# methods - graph-based gene module detection ----------------------------------

#' Identify correlation-based gene modules via graphs.
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
#' @importFrom magrittr `%$%`
#' @importFrom zeallot `%<-%`
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

  graph <- if(S7::prop(bulk_coexp, "params")[["detection_method"]] == "correlation-based") {
    get_cor_graph_single(
      bulk_coexp,
      kernel_bandwidth = kernel_bandwidth,
      min_affinity = min_affinity,
      .verbose = .verbose
    )
  }

  resolutions <- exp(seq(
    log(min_res),
    log(max_res),
    length.out = number_res
  ))

  if (.verbose) message(sprintf("Iterating through %i resolutions", length(resolutions)))

  if(parallel) {
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
    .progress = .verbose,
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


