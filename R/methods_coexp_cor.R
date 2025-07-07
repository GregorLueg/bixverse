# methods - simple correlations ------------------------------------------------

#' @title Prepare correlation-based module detection
#'
#' @description
#' This function will calculate the correlation coefficients between the genes,
#' using the highly variable genes (if available, otherwise the function will
#' use the raw data). The data will be stored in a memory-efficient format
#' in the properties of the class.
#'
#' @param object The class, see [bixverse::bulk_coexp()]. Ideally, you
#' should run [bixverse::preprocess_bulk_coexp()] before applying this function.
#' @param cor_method String. Option of `c("pearson", "spearman")`.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @return The class with added data to the properties for subsequent usage.
#'
#' @export
cor_module_processing <- S7::new_generic(
  name = "cor_module_processing",
  dispatch_args = "object",
  fun = function(
    object,
    cor_method = c("pearson", "spearman"),
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)


#' @export
#'
#' @method cor_module_processing bulk_coexp
S7::method(cor_module_processing, bulk_coexp) <- function(
  object,
  cor_method = c("pearson", "spearman"),
  .verbose = TRUE
) {
  # Checks
  checkmate::assertClass(object, "bixverse::bulk_coexp")
  checkmate::assertChoice(cor_method, c("pearson", "spearman"))
  checkmate::qassert(.verbose, "B1")

  # Function body
  if (purrr::is_empty(S7::prop(object, "processed_data")[["processed_data"]])) {
    warning("No pre-processed data found. Defaulting to the raw data.")
    target_mat <- S7::prop(object, "raw_data")
  } else {
    target_mat <- S7::prop(object, "processed_data")[["processed_data"]]
  }

  spearman <- if (cor_method == "pearson") {
    if (.verbose) {
      message("Using Pearson correlations.")
    }
    FALSE
  } else {
    if (.verbose) {
      message("Using Spearman correlations.")
    }
    TRUE
  }

  # Calculate the upper triangle of correlation matrix
  cor_diagonal <- rs_cor_upper_triangle(
    target_mat,
    spearman = spearman,
    shift = 1L
  )

  # Save data to memory friendly R6 class
  cor_data <- upper_triangular_sym_mat$new(
    values = cor_diagonal,
    features = colnames(target_mat),
    shift = 1L
  )

  correlation_params <- list(spearman = spearman, type = "simple")

  S7::prop(object, "processed_data")[["correlation_res"]] <- cor_data
  S7::prop(object, "params")[["correlation_params"]] <- correlation_params
  S7::prop(object, "params")["detection_method"] <- "correlation-based"

  return(object)
}

# methods - TOM ----------------------------------------------------------------

#' @title Update the correlation matrix to a TOM
#'
#' @description
#' This function will update the correlation matrix to a topological overlap
#' matrix. It defaults to `"v2"` and the signed version, please see
#' [bixverse::calculate_tom()] for details.
#'
#' @param object The class, see [bixverse::bulk_coexp()]. You need to have
#' applied [bixverse::cor_module_processing()] before applying this function.
#' @param signed Boolean. Do you want to use the signed or unsigned version.
#' Defaults to `TRUE`.
#' @param version String. One of `c("v2", "v1")`. Defaults to `"v2"`.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @return The class with added data to the properties for subsequent usage.
#'
#' @export
cor_module_tom <- S7::new_generic(
  name = "cor_module_tom",
  dispatch_args = "object",
  fun = function(
    object,
    signed = TRUE,
    version = c("v2", "v1"),
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @export
#'
#' @method cor_module_tom bulk_coexp
S7::method(cor_module_tom, bulk_coexp) <- function(
  object,
  signed = TRUE,
  version = c("v2", "v1"),
  .verbose = TRUE
) {
  version <- match.arg(version)

  # checks
  checkmate::assertClass(object, "bixverse::bulk_coexp")
  checkmate::qassert(signed, "B1")
  checkmate::assertChoice(version, c("v2", "v1"))
  checkmate::qassert(.verbose, "B1")

  # early return
  detection_method <- S7::prop(object, "params")[["detection_method"]]
  if (
    is.null(detection_method) ||
      detection_method != "correlation-based"
  ) {
    warning(
      paste(
        "This class does not seem to be set for correlation-based module detection.",
        "Returning class as is."
      )
    )
    return(object)
  }

  # pull out the correlation results
  cor_res <- S7::prop(object, "processed_data")$correlation_res
  cor_mat <- cor_res$get_sym_matrix(.verbose = .verbose)

  if (.verbose) {
    message("Replacing the correlation matrix with a TOM.")
  }

  features <- rownames(cor_mat)
  tom_mat <- rs_tom(x = cor_mat, tom_type = version, signed = signed)
  tom_vec <- rs_dense_to_upper_triangle(tom_mat, 1L)

  tom_res <- upper_triangular_sym_mat$new(
    values = tom_vec,
    features = features,
    shift = 1L
  )

  S7::prop(object, "processed_data")[["correlation_res"]] <- tom_res
  S7::prop(object, "params")[["correlation_params"]][["TOM"]] <- TRUE

  return(object)
}

# methods - differential correlations ------------------------------------------

#' @title Prepare differential correlation-based module detection
#'
#' @description
#' This function will calculate the differential correlation between the stored
#' data set in the class and another background data set. To do so, it uses a
#' Fisher transformation of the correlation coefficients and calculates a Z
#' score based on the delta. The function will automatically subset into shared
#' features between the two data sets.
#'
#' @param object The class, see [bixverse::bulk_coexp()]. Ideally, you
#' should run [bixverse::preprocess_bulk_coexp()] before applying this function.
#' @param background_mat Numerical matrix. The background data set.
#' @param cor_method String. Option of `c("pearson", "spearman")`.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @return The class with added data to the properties for subsequent usage.
#'
#' @export
diffcor_module_processing <- S7::new_generic(
  name = "diffcor_module_processing",
  dispatch_args = "object",
  fun = function(
    object,
    background_mat,
    cor_method = c("pearson", "spearman"),
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @export
#'
#' @method diffcor_module_processing bulk_coexp
S7::method(diffcor_module_processing, bulk_coexp) <- function(
  object,
  background_mat,
  cor_method = c("pearson", "spearman"),
  .verbose = TRUE
) {
  # Checks
  checkmate::assertClass(object, "bixverse::bulk_coexp")
  checkmate::assertMatrix(background_mat, mode = "numeric")
  checkmate::assertChoice(cor_method, c("pearson", "spearman"))
  checkmate::qassert(.verbose, "B1")

  # Function
  if (purrr::is_empty(S7::prop(object, "processed_data")[["processed_data"]])) {
    warning("No pre-processed data found. Defaulting to the raw data.")
    target_mat <- S7::prop(object, "raw_data")
  } else {
    target_mat <- S7::prop(object, "processed_data")[["processed_data"]]
  }

  spearman <- if (cor_method == "pearson") {
    if (.verbose) {
      message("Using Pearson correlations.")
    }
    FALSE
  } else {
    if (.verbose) {
      message("Using Spearman correlations.")
    }
    TRUE
  }

  features <- colnames(target_mat)
  shared_features <- intersect(colnames(target_mat), colnames(background_mat))

  # Early return if there are no shared features
  if (length(shared_features) == 0) {
    warning("No shared features identified. Returning class as is.")
    return(object)
  }
  if (.verbose) {
    message(
      sprintf(
        "A total of %i shared features were identified and used for differential correlation.",
        length(shared_features)
      )
    )
  }

  target_mat <- target_mat[, shared_features]
  background_mat <- background_mat[, shared_features]

  combined_mad_df <- list(
    feature = shared_features,
    MAD_target = matrixStats::colMads(target_mat),
    MAD_background = matrixStats::colMads(background_mat)
  ) %>%
    data.table::setDT()

  diff_cor_res <- rs_differential_cor(
    target_mat,
    background_mat,
    spearman = spearman
  )

  cor_data <- upper_triangle_diffcor_mat$new(
    diff_cor_res = diff_cor_res,
    features = shared_features
  )

  correlation_params <- list(
    spearman = spearman,
    type = "differential correlation",
    no_intersecting_features = length(shared_features)
  )

  S7::prop(object, "processed_data")[["differential_cor_res"]] <- cor_data
  S7::prop(object, "processed_data")[[
    "differential_cor_feature_meta"
  ]] <- combined_mad_df
  S7::prop(object, "params")[["correlation_params"]] <- correlation_params
  S7::prop(object, "params")[
    "detection_method"
  ] <- "differential correlation-based"

  object
}


# methods - graph-based gene module detection ----------------------------------

#' @title Iterate through different epsilon parameters
#'
#' @description
#' This functions iterates through a set of provided epsilons and checks for
#' each one to what extend the underlying affinity matrix will follow a power
#' law distribution.
#'
#' @param object The class, see [bixverse::bulk_coexp()]. You need to run
#' [bixverse::cor_module_processing()] before running this function.
#' @param rbf_func The type of RBF function you want to apply. A choice of
#' `c('bump', 'gaussian', 'inverse_quadratic')`.
#' @param epsilons Vector of floats. The different epsilon parameters you
#' would like to run.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @return The class with added data to the properties for subsequent usage.
#'
#' @export
cor_module_check_epsilon <- S7::new_generic(
  name = "cor_module_check_epsilon",
  dispatch_args = "object",
  fun = function(
    object,
    rbf_func = c('bump', 'gaussian', 'inverse_quadratic'),
    epsilons = c(0.25, seq(from = 0.5, to = 5, by = 0.5)),
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @export
#'
#' @importFrom magrittr `%>%`
#' @importFrom zeallot `%->%`
#' @import data.table
#'
#' @method cor_module_check_epsilon bulk_coexp
S7::method(cor_module_check_epsilon, bulk_coexp) <- function(
  object,
  rbf_func = c('bump', 'gaussian', 'inverse_quadratic'),
  epsilons = c(0.25, seq(from = 0.5, to = 5, by = 0.5)),
  .verbose = TRUE
) {
  # Checks
  checkmate::assertClass(object, "bixverse::bulk_coexp")
  checkmate::qassert(epsilons, "R+")
  checkmate::qassert(.verbose, "B1")
  checkmate::assertChoice(rbf_func, c('bump', 'gaussian', 'inverse_quadratic'))

  detection_method <- S7::prop(object, "params")[["detection_method"]]

  # Early return
  if (
    is.null(detection_method) ||
      detection_method != "correlation-based"
  ) {
    warning(
      paste(
        "This class does not seem to be set for correlation-based module detection.",
        "Returning class as is."
      )
    )
    return(object)
  }

  # Pull out the correlation results
  cor_res <- S7::prop(object, "processed_data")$correlation_res
  c(cor_vector, features, n_features, shift) %<-% cor_res$get_data()

  # Prepare everything for iterating through the epsilons
  epsilons <- sort(epsilons, decreasing = TRUE)
  dist_vec <- 1 - abs(cor_vector)
  dist_vec <- data.table::fifelse(dist_vec < 0, 0, dist_vec)

  if (.verbose) {
    message(sprintf("Testing %i epsilons.", length(epsilons)))
  }

  epsilon_data <- rs_rbf_iterate_epsilons(
    dist = dist_vec,
    epsilon_vec = epsilons,
    original_dim = n_features,
    shift = shift,
    rbf_type = rbf_func
  )

  r_square_vals <- apply(epsilon_data, 2, scale_free_fit)

  r_square_data <- list(epsilon = epsilons, r2_vals = r_square_vals) %>%
    data.table::setDT()

  S7::prop(object, "outputs")[["epsilon_data"]] <- r_square_data

  return(object)
}


#' @title Iterate through Leiden resolutions for graph-based community detection.
#'
#' @description
#' This function will identify gene modules based on affinity graphs from the
#' single correlation or differential correlation methods. Briefly, in the case
#' of single correlation, the graph is generated based on the absolute
#' correlation coefficients that are subjected to a Gaussian affinity kernel.
#' This reduces spurious correlations and leaves a sparsely connected graph.
#' In the case of differential correlations, the graph is generated based on
#' significant differential correlations if one of the two correlations reached
#' the defined minimum thresholds.\cr
#' Subsequently, Leiden community detection is applied on the respective graph
#' through a range of resolutions that the user can define. The function then
#' returns meta information about the resolutions (which can also be plotted) to
#' identify the best suitable resolution parameter to identify co-expression modules.
#'
#' @param object The class, see [bixverse::bulk_coexp()].
#' @param resolution_params List. Parameters for the resolution search, see
#' [bixverse::params_graph_resolution()]. Contains:
#' \itemize{
#'  \item min_res - Float. Minimum resolution to test.
#'  \item max_res - Float. Maximum resolution to test.
#'  \item number_res - Integer. Number of resolutions to test between the
#'  `max_res` and `min_res.`
#' }
#' @param graph_params List. Parameters for the generation of the (differential)
#' correlation graph, see [bixverse::params_cor_graph()]. Contains:
#' \itemize{
#'  \item Epsilon - Defines the epsilon parameter for the radial basis
#'  function. Defaults to 1, but should be ideally optimised.
#'  \item min_cor - Float. Minimum absolute correlation that needs to be
#'  observed in either data set. Only relevant for differential correlation-based
#'  graphs.
#'  \item fdr_threshold - Float. Maximum FDR for the differential correlation
#'  p-value.
#'  \item verbose - Boolean. Controls verbosity of the graph generation.
#' }
#' @param random_seed Integer. Random seed.
#' @param min_genes Integer. Minimum number of genes that should be in a
#' community.
#' @param parallel Boolean. Parallelise the Leiden clustering.
#' @param max_workers Optional Integer. Number of cores to use if parallel is
#' set to `TRUE`. If set to `NULL` it will automatically detect the number
#' of cores.
#' @param .verbose Controls the verbosity of the function.
#'
#' @return The class with added data to the properties.
#'
#' @export
cor_module_graph_check_res <- S7::new_generic(
  name = "cor_module_graph_check_res",
  dispatch_args = "object",
  fun = function(
    object,
    resolution_params = params_graph_resolution(),
    graph_params = params_cor_graph(),
    random_seed = 123L,
    min_genes = 10L,
    parallel = TRUE,
    max_workers = NULL,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)


#' @export
#'
#' @importFrom magrittr `%>%`
#' @importFrom zeallot `%->%`
#' @importFrom future plan multisession sequential
#' @import data.table
#'
#' @method cor_module_graph_check_res bulk_coexp
S7::method(cor_module_graph_check_res, bulk_coexp) <- function(
  object,
  resolution_params = params_graph_resolution(),
  graph_params = params_cor_graph(),
  random_seed = 123L,
  min_genes = 10L,
  parallel = TRUE,
  max_workers = NULL,
  .verbose = TRUE
) {
  combined_id <- . <- good_clusters <- N <- graph <- NULL

  # Checks
  checkmate::assertClass(object, "bixverse::bulk_coexp")
  assertCorGraphParams(graph_params)
  assertGraphResParams(resolution_params)
  checkmate::qassert(min_genes, "I1")
  checkmate::qassert(parallel, "B1")
  checkmate::qassert(max_workers, c("I1", "0"))
  checkmate::qassert(.verbose, "B1")

  detection_method <- S7::prop(object, "params")[["detection_method"]]

  # Early return
  if (
    is.null(detection_method) ||
      !detection_method %in%
        c("correlation-based", "differential correlation-based")
  ) {
    warning(
      paste(
        "This class does not seem to be set for correlation-based module detection",
        "Returning class as is."
      )
    )
    return(object)
  }

  c(graph, graph_params) %<-%
    with(
      graph_params,
      switch(
        detection_method,
        "correlation-based" = get_cor_graph(
          object = object,
          epsilon = epsilon,
          .verbose = verbose
        ),
        "differential correlation-based" = get_diffcor_graph(
          object = object,
          min_cor = min_cor,
          fdr_threshold = fdr_threshold,
          .verbose = verbose
        )
      )
    )

  resolutions <- with(
    resolution_params,
    exp(seq(log(min_res), log(max_res), length.out = number_res))
  )

  if (.verbose) {
    message(sprintf("Iterating through %i resolutions", length(resolutions)))
  }

  if (parallel) {
    if (is.null(max_workers)) {
      max_workers <- get_cores()
    }
    if (.verbose) {
      message(sprintf("Using parallel computation over %i cores.", max_workers))
    }

    # future plan funkiness
    assign(".temp_workers", max_workers, envir = .GlobalEnv)
    on.exit(rm(".temp_workers", envir = .GlobalEnv))

    plan(future::multisession(workers = .temp_workers))
  } else {
    if (.verbose) {
      message("Using sequential computation.")
    }
    future::plan(future::sequential())
  }

  community_df_res <- furrr::future_map(
    resolutions,
    \(res) {
      set.seed(random_seed)
      community <- igraph::cluster_leiden(
        graph,
        objective_function = "modularity",
        resolution = res,
        n_iterations = 5L
      )

      modularity <- igraph::modularity(
        x = graph,
        membership = community$membership
      )

      community_df <- data.table::data.table(
        resolution = res,
        node_name = community$names,
        membership = community$membership,
        modularity = modularity
      )
    },
    .progress = .verbose,
    .options = furrr::furrr_options(seed = TRUE)
  ) %>%
    data.table::rbindlist(.)

  # To make the message trace prettier, if set to verbose
  if (.verbose) {
    message("")
  }

  future::plan(future::sequential())
  gc()

  community_df_res[, combined_id := sprintf("id_%s_%s", resolution, membership)]

  cluster_summary <- community_df_res[, .N, combined_id] %>%
    .[, good_clusters := N >= min_genes] %>%
    data.table::merge.data.table(
      .,
      unique(community_df_res[, c("resolution", "combined_id")]),
      by.x = "combined_id",
      by.y = "combined_id",
    ) %>%
    .[,
      .(
        good_clusters = sum(good_clusters),
        avg_size = mean(N),
        max_size = max(N)
      ),
      resolution
    ]

  resolution_results <- community_df_res[,
    .(
      no_clusters = length(unique(membership)),
      modularity = unique(modularity)
    ),
    resolution
  ] %>%
    data.table::merge.data.table(
      .,
      cluster_summary,
      by.x = "resolution",
      by.y = "resolution"
    )

  # Assign stuff
  S7::prop(object, "params")[["correlation_graph"]] <- graph_params
  S7::prop(object, "outputs")[["resolution_results"]] <- resolution_results
  S7::prop(object, "outputs")[["cor_graph"]] <- graph

  return(object)
}

#' @title Identify correlation-based gene modules via graphs
#'
#' @description
#' This function leverages graph-based clustering to identify gene co-expression
#' modules. The class has the option to sub-cluster large communities within
#' their respective sub graphs, akin to the approach taken by Barrio-Hernandez,
#' et al.
#'
#' @param object The class, see [bixverse::bulk_coexp()].
#' @param resolution The Leiden resolution parameter you wish to use. If NULL,
#' it will use the optimal one identified by [bixverse::cor_module_graph_check_res()].
#' If nothing can be found, will default to 1.
#' @param min_size Integer. Minimum size of the communities.
#' @param max_size Integer. Maximum size of the communities.
#' @param subclustering Boolean. Shall after a first clustering communities that
#' are too large be further sub clustered. Defaults to `TRUE`.
#' @param random_seed Integer. Random seed.
#' @param .graph_params List. Parameters for the generation of the (differential)
#' correlation graph, see [bixverse::params_cor_graph()]. Contains:
#' \itemize{
#'  \item Epsilon - Defines the epsilon parameter for the radial basis
#'  function. Defaults to 2, but should be ideally optimised.
#'  \item min_cor - Float. Minimum absolute correlation that needs to be
#'  observed in either data set. Only relevant for differential correlation-based
#'  graphs.
#'  \item fdr_threshold - Float. Maximum FDR for the differential correlation
#'  p-value.
#'  \item verbose - Boolean. Controls verbosity of the graph generation.
#' }
#' This parameter is only relevant if you did *not* run
#' [bixverse::cor_module_graph_check_res()].
#' @param .max_iters Integer. If sub clustering is set to `TRUE`, what shall be the
#' maximum number of iterations. Defaults to 100L.
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @return The class with added data to the properties.
#'
#' @references Barrio-Hernandez, et al., Nat Genet, 2023.
#'
#' @export
cor_module_graph_final_modules <- S7::new_generic(
  name = "cor_module_graph_final_modules",
  dispatch_args = "object",
  fun = function(
    object,
    resolution = NULL,
    min_size = 10L,
    max_size = 500L,
    subclustering = TRUE,
    random_seed = 123L,
    .graph_params = params_cor_graph(),
    .max_iters = 100L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @export
#'
#' @importFrom magrittr `%>%`
#' @importFrom zeallot `%->%`
#' @import data.table
#'
#' @method cor_module_graph_final_modules bulk_coexp
S7::method(cor_module_graph_final_modules, bulk_coexp) <- function(
  object,
  resolution = NULL,
  min_size = 10L,
  max_size = 500L,
  subclustering = TRUE,
  random_seed = 123L,
  .graph_params = params_cor_graph(),
  .max_iters = 100L,
  .verbose = TRUE
) {
  graph_params <- graph <- modularity <- . <- N <- cluster_id <- node_id <- NULL

  # Checks
  checkmate::assertClass(object, "bixverse::bulk_coexp")
  checkmate::qassert(resolution, c("0", "N1"))
  checkmate::qassert(min_size, "I1")
  checkmate::qassert(max_size, "I1")
  checkmate::qassert(subclustering, "B1")
  checkmate::qassert(random_seed, "I1")
  assertCorGraphParams(.graph_params)
  checkmate::qassert(.max_iters, "I1")
  checkmate::qassert(.verbose, "B1")

  detection_method <- S7::prop(object, "params")[["detection_method"]]

  # Early return
  if (
    is.null(detection_method) &&
      detection_method %in%
        c("correlation-based", "differential correlation-based")
  ) {
    warning(
      paste(
        "This class does not seem to be set for correlation-based module detection",
        "Returning class as is."
      )
    )
    return(object)
  }

  # Get the graph
  if (is.null(S7::prop(object, "outputs")[["cor_graph"]])) {
    # Deal with the case a graph was not yet generated...
    warning(
      paste(
        "No correlation graph found. Did you run cor_module_graph_check_res()?",
        "Generating correlation graph based on standard parameters.",
        sep = "\n"
      )
    )

    c(graph, graph_params) %<-%
      with(
        .graph_params,
        switch(
          detection_method,
          "correlation-based" = get_cor_graph(
            object = object,
            kernel_bandwidth = kernel_bandwidth,
            min_affinity = min_affinity,
            .verbose = verbose
          ),
          "differential correlation-based" = get_diffcor_graph(
            object = object,
            min_cor = min_cor,
            fdr_threshold = fdr_threshold,
            .verbose = verbose
          )
        )
      )

    S7::prop(object, "params")[["correlation_graph"]] <- graph_params
    S7::prop(object, "outputs")[["cor_graph"]] <- graph
  } else {
    graph <- S7::prop(object, "outputs")[["cor_graph"]]
  }

  # Final resolution
  if (is.null(resolution)) {
    resolution_results <- S7::prop(object, "outputs")[["resolution_results"]]
    final_resolution <- if (!is.null(resolution_results)) {
      if (.verbose) {
        message("Using resolution with best modularity.")
      }
      resolution_results[modularity == max(modularity), resolution]
    } else {
      warning(
        paste(
          "No resolution results found and none provided.",
          "Will default to a resolution of 1."
        )
      )
      1
    }
  }

  # Do a first clustering
  set.seed(random_seed)

  final_gene_communities <- igraph::cluster_leiden(
    graph,
    objective_function = "modularity",
    resolution = final_resolution,
    n_iterations = 5L
  )

  clusters_df <- data.table::data.table(
    node_id = final_gene_communities$names,
    cluster_id = final_gene_communities$membership
  )

  node_frequency <- clusters_df[, .N, .(cluster_id)]

  # If sub clustering is active, do that
  if (subclustering) {
    if (.verbose) {
      message(
        "Sub-clustering larger communities until they are below max_size."
      )
    }
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

        sub_graph_l <- igraph::subgraph(
          graph,
          data.table::chmatch(nodes_in_cluster, igraph::V(graph)$name)
        )

        # Restarting at a very small resolution
        clusters_red <- igraph::cluster_leiden(
          sub_graph_l,
          resolution = 0.1 + l * 0.05
        )

        subclusters <- data.table(
          node_id = clusters_red$names,
          cluster_id = clusters_red$membership
        )

        subclusters_frequency <- subclusters[, .N, .(cluster_id)]
        clusters_small_enough <- subclusters_frequency[
          N <= max_size,
          cluster_id
        ]

        good_clusters <- subclusters[cluster_id %in% clusters_small_enough] %>%
          dplyr::mutate(
            cluster_id = paste0(
              i,
              paste(rep("sub", l), collapse = ""),
              cluster_id
            )
          )

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
  final_communities <- node_frequency_updated[
    N <= max_size &
      N >= min_size,
    cluster_id
  ]
  final_clusters_filtered <- final_clusters[cluster_id %in% final_communities]

  cluster_name_prettifier <- setNames(
    paste(
      "cluster",
      seq_along(
        unique(final_clusters_filtered$cluster_id)
      ),
      sep = "_"
    ),
    unique(final_clusters_filtered$cluster_id)
  )

  final_clusters_filtered[, cluster_id := cluster_name_prettifier[cluster_id]]

  results_param <- list(
    resolution = resolution,
    seed = random_seed,
    min_size = min_size,
    max_size = max_size,
    max_iters = .max_iters
  )

  S7::prop(object, "final_results") <- final_clusters_filtered
  S7::prop(object, "params")[["module_final_gen"]] <- results_param

  return(object)
}

# methods - coremo -------------------------------------------------------------

#' @title Generates CoReMo-based gene modules
#'
#' @description
#' This function creates gene modules, based on the framework from Srivastava
#' et al., 2018. Briefly, it applies an RBF function to the correlation matrix
#' to reduce the impact of weak correlations, leverages hierarchical clustering
#' for clustering the genes. It optimises the R2 (cor^2) within the clusters to
#' identify the optimal cut. Gene modules with low R2 are being considered as
#' the 'junk module'.
#'
#' @param object The class, see [bixverse::bulk_coexp()].
#' @param coremo_params List. Parameters for the generation of the CoReMo
#' modules, see [bixverse::params_coremo()]. Contains:
#' \itemize{
#'  \item epsilon - Float. The epsilon parameter for the RBF. You can optimise
#'  that one with [bixverse::cor_module_check_epsilon()]. Defaults to `2`.
#'  \item k_min - Integer. Minimum number of cuts. Defaults to `2L`.
#'  \item k_max - Integer. Maximum number of cuts. Defaults to `150L`.
#'  \item min_size - Optional integer. Minimum size of the clusters. If
#'  provided, smaller clusters will be merged by eigengene similarity.
#'  \item junk_module_threshold - Float. Minimum R2 median  value for a module
#'  to not be considered a junk module. Defaults to `0.05`.
#'  \item rbf_func - String. Type of RBF function to apply. Defaults to
#'  `"gaussian"`.
#'  \item cor_method - String. Type of correlation method to use for merging
#'  the smaller cluster. Defaults to `"spearman"`.
#' }
#' @param seed Integer. Random seed for reproducibility purposes.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @return The class with added data to the properties for subsequent usage.
#'
#' @references Srivastava, et al., Nat. Commun., 2018
#'
#' @export
cor_module_coremo_clustering <- S7::new_generic(
  name = "cor_module_coremo_clustering",
  dispatch_args = "object",
  fun = function(
    object,
    coremo_params = params_coremo(),
    seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @export
#'
#' @importFrom magrittr `%>%`
#' @importFrom zeallot `%<-%`
#' @import data.table
#'
#' @method cor_module_coremo_clustering bulk_coexp
S7::method(cor_module_coremo_clustering, bulk_coexp) <- function(
  object,
  coremo_params = params_coremo(),
  seed = 42L,
  .verbose = TRUE
) {
  # Out of scope
  gradient_change <- r2med <- cluster_id <- . <- NULL
  # Checks
  checkmate::assertClass(object, "bixverse::bulk_coexp")
  assertCoReMoParams(coremo_params)
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, "B1")

  detection_method <- S7::prop(object, "params")[["detection_method"]]

  # Early return
  if (
    is.null(detection_method) ||
      detection_method != "correlation-based"
  ) {
    warning(
      paste(
        "This class does not seem to be set for correlation-based module detection.",
        "Returning class as is."
      )
    )
    return(object)
  }

  cor_res <- S7::prop(object, "processed_data")[["correlation_res"]]
  cor_mat <- cor_res$get_sym_matrix(.verbose = .verbose)

  aff_mat <- with(
    coremo_params,
    rs_rbf_function_mat(
      x = 1 - abs(cor_mat),
      epsilon = epsilon,
      rbf_type = rbf_func
    )
  )
  dist_mat <- 1 - aff_mat

  if (.verbose) {
    message("Generating the hierarchical clustering.")
  }
  tree <- fastcluster::hclust(as.dist(dist_mat), method = "ward.D")

  if (.verbose) {
    message("Identifying optimal number of cuts.")
  }
  optimal_cuts <- with(
    coremo_params,
    tree_cut_iter(
      tree = tree,
      cor_mat = cor_mat,
      dist_mat = dist_mat,
      k_min = k_min,
      k_max = k_max,
      min_size = min_size,
      cor_method = cor_method,
      seed = seed
    )
  )

  c(inflection_idx, gradient_change) %<-%
    get_inflection_point(
      optimal_cuts$k,
      optimal_cuts$R2_weighted_median,
      span = 0.25
    )

  optimal_cuts[, gradient_change := c(0, gradient_change)]

  if (.verbose) {
    message("Finalising CoReMo clusters.")
  }
  final_clusters <- with(
    coremo_params,
    coremo_tree_cut(
      tree = tree,
      k = as.integer(inflection_idx),
      dist_mat = dist_mat,
      cor_method = cor_method
    )
  ) %>%
    `names<-`(rownames(cor_mat))

  cluster_list <- split(
    names(final_clusters),
    final_clusters
  )

  final_quality <- rs_coremo_quality(
    cluster_genes = cluster_list,
    cor_mat = cor_mat,
    row_names = rownames(cor_mat),
    seed = seed
  ) %>%
    data.table::setDT() %>%
    .[, cluster_id := names(cluster_list)]

  junk_modules <- final_quality[
    r2med <= coremo_params$junk_module_threshold,
    cluster_id
  ]

  module_dt <- data.table::as.data.table(
    stack(cluster_list)
  ) %>%
    data.table::setnames(
      old = c("values", "ind"),
      new = c("gene", "cluster_id")
    ) %>%
    .[!cluster_id %in% junk_modules] %>%
    merge(., final_quality, by = "cluster_id")

  coremo_param <- with(
    coremo_params,
    list(
      epsilon = epsilon,
      k_min = k_min,
      k_max = k_max,
      cor_method = cor_method,
      min_size = min_size,
      seed = seed,
      rbf_func = rbf_func,
      junk_module_threshold = junk_module_threshold,
      inflection_idx = inflection_idx
    )
  )

  # Assign the objects
  S7::prop(object, "params")[["coremo"]] <- coremo_param
  S7::prop(object, "outputs")[["optimal_cuts"]] <- optimal_cuts
  S7::prop(object, "outputs")[["final_modules"]] <- module_dt
  S7::prop(object, "outputs")[["tree"]] <- tree

  return(object)
}

#' @title Assesses CoReMo-based gene module stability
#'
#' @description
#' The function assesses the stability of the CoReMo modules, leveraging a
#' leave-one-out sample method. In each iteration, one of the samples is left
#' out and the clustering is redone. Subsequently, the Jaccard similarity (as
#' a surrogate for membership stability) is calculated for each feature across
#' the different resamplings and added to the `final_modules` data.table.
#'
#' @param object The class, see [bixverse::bulk_coexp()].
#' @param chunk_size Integer. Chunk size in which to process the data. Defaults
#' to `15L`, i.e., 15 samples are being processed in one go. You can use bigger
#' values here, but be aware of memory pressure.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @return The class with added data to the properties for subsequent usage.
#'
#' @references Srivastava, et al., Nat. Commun., 2018; Francois, Romagnolo,
#' et al., Nat. Commun., 2024.
#'
#' @export
cor_module_coremo_stability <- S7::new_generic(
  name = "cor_module_coremo_stability",
  dispatch_args = "object",
  fun = function(
    object,
    chunk_size = 15L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @method cor_module_coremo_stability bulk_coexp
S7::method(cor_module_coremo_stability, bulk_coexp) <- function(
  object,
  chunk_size = 15L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertClass(object, "bixverse::bulk_coexp")
  checkmate::qassert(chunk_size, "I1")
  checkmate::qassert(.verbose, "B1")

  # function
  # early return
  if (is.null(S7::prop(object, "params")[["coremo"]])) {
    warning(paste(
      "No CoReMo-data found. Did you run cor_module_coremo_clustering()?",
      "Returning class as is."
    ))
    return(object)
  }
  if (purrr::is_empty(S7::prop(object, "processed_data")[["processed_data"]])) {
    warning("No pre-processed data found. Defaulting to the raw data.")
    data_mat <- S7::prop(object, "raw_data")
  } else {
    data_mat <- S7::prop(object, "processed_data")[["processed_data"]]
  }

  # pull out the needed parameters
  coremo_params <- S7::prop(object, "params")[["coremo"]]
  final_modules <- S7::prop(object, "outputs")[["final_modules"]]

  total_samples <- seq_len(nrow(data_mat))
  no_total_samples <- length(total_samples)
  groups <- ceiling(seq_along(total_samples) / chunk_size)
  chunks <- split(total_samples, groups)

  if (.verbose) {
    message(sprintf(
      "Running the leave-one-out stability assessment over %i chunks.",
      length(chunks)
    ))
  }

  all_results <- vector(mode = "list", length = length(chunks))

  for (i in seq_along(chunks)) {
    indices <- chunks[[i]]
    chunk_res <- rs_coremo_stability(
      data = data_mat,
      indices = indices,
      epsilon = coremo_params$epsilon,
      rbf_type = coremo_params$rbf_func,
      spearman = TRUE
    )

    leave_one_out_clustering <- purrr::map(chunk_res, \(chunk) {
      dist_obj <- create_dist_obj(
        x = chunk,
        size = ncol(data_mat)
      )
      new_tree <- fastcluster::hclust(dist_obj, method = "ward.D")
      clusters <- cutree(new_tree, k = coremo_params$inflection_idx)

      clusters
    })

    all_results[[i]] <- leave_one_out_clustering

    message_txt <- sprintf(
      "Chunk %i out of %i: Processed %i samples out of a total of %i samples.",
      i,
      length(chunks),
      ifelse(
        chunk_size * i < no_total_samples,
        chunk_size * i,
        no_total_samples
      ),
      no_total_samples
    )

    if (.verbose) message(message_txt)
  }

  all_results <- purrr::flatten(all_results)

  cluster_mat <- do.call(cbind, all_results) %>%
    `rownames<-`(colnames(data_mat))

  if (.verbose) {
    message(
      "Assessing stability of gene membership within leave-one-out resamling."
    )
  }

  stability <- rs_cluster_stability(cluster_mat[final_modules$gene, ])

  final_modules[, c("stability", "std_stability") := stability]

  S7::prop(object, "outputs")[["final_modules"]] <- final_modules

  return(object)
}

# methods - helpers ------------------------------------------------------------

## power law calculations ------------------------------------------------------

#' Calculate the goodness of fit for a power law distribution.
#'
#' @param k Numeric vector. The vector of node degrees.
#' @param breaks Integer. Number of breaks for fitting the data.
#' @param plot Boolean. Shall the log-log plot be generated.
#'
#' @returns The R2 value of of the goodness of fit.
scale_free_fit <- function(k, breaks = 50L, plot = FALSE) {
  # Visible global function stuff...
  lm <- NULL
  # Checks
  checkmate::qassert(k, "R>=50")
  checkmate::qassert(breaks, "I1")
  checkmate::qassert(plot, "B1")
  # Deal with stupid case that someone supplies something small here
  if (breaks > length(k)) {
    breaks <- ceiling(k / 10L)
  }
  k_discrete <- cut(k, breaks)
  dk <- tapply(k, k_discrete, mean)
  dk_p <- tapply(k, k_discrete, length) / length(k)
  log_dk <- log10(dk)
  log_dk_p <- log10(dk_p)

  if (plot) {
    plot(x = log_dk, y = log_dk_p)
  }

  summary(lm(log_dk ~ log_dk_p))$r.squared
}

## graph generation ------------------------------------------------------------

#' @title Get correlation-based graph
#'
#' @description
#' Helper function to get a correlation-based igraph from the class
#'
#' @param object The class, see [bixverse::bulk_coexp()].
#' @param epsilon Float. The epsilon parameter for the RBF function, in this
#' case the bump function.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @return A list with the following elements:
#' \itemize{
#'  \item graph - The igraph
#'  \item params - A list that contains the parameters of the graph generation
#'  and general graph information (node, edge numbers).
#' }
#'
#' @export
get_cor_graph <- S7::new_generic(
  name = "get_cor_graph",
  dispatch_args = "object",
  fun = function(object, epsilon, .verbose) {
    S7::S7_dispatch()
  }
)


#' @export
#' @method get_cor_graph bulk_coexp
S7::method(get_cor_graph, bulk_coexp) <- function(object, epsilon, .verbose) {
  # Scope checks...
  . <- cor_abs <- affinity <- NULL
  # Checks
  checkmate::assertClass(object, "bixverse::bulk_coexp")
  checkmate::qassert(epsilon, "R1")
  checkmate::qassert(.verbose, "B1")
  # Function body
  cor_res <- S7::prop(object, "processed_data")$correlation_res
  graph_df <- cor_res$get_data_table(.verbose = .verbose) %>%
    .[, cor_abs := abs(cor)] %>%
    .[, dist := 1 - cor_abs] %>%
    .[, dist := data.table::fifelse(dist < 0, 0, dist)] %>%
    .[,
      affinity := rs_rbf_function(
        x = dist,
        epsilon = epsilon,
        rbf_type = "bump"
      )
    ] %>%
    .[affinity > 0] %>%
    .[, c("feature_a", "feature_b", "affinity")] %>%
    data.table::setnames(
      .,
      old = c("feature_a", "feature_b", "affinity"),
      new = c("from", "to", "weight")
    )

  if (.verbose) {
    message(sprintf(
      "Generating correlation-based graph with %i edges.",
      nrow(graph_df)
    ))
  }

  graph <- igraph::graph_from_data_frame(graph_df, directed = FALSE)
  graph <- igraph::simplify(graph)

  graph_params <- list(
    epsilon = epsilon,
    no_nodes = length(igraph::V(graph)),
    no_edges = length(igraph::E(graph))
  )

  list(graph = graph, params = graph_params)
}

#' @title Get differential correlation-based graph
#'
#' @description
#' Helper function to get a differential correlation-based igraph from the class
#'
#' @param object The class, see [bixverse::bulk_coexp()].
#' @param min_cor Float. The minimum absolute correlation that needs to be
#' present in either data set.
#' @param fdr_threshold Float. The maximum FDR that is tolerated for the
#' generation of the graph.
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @return A list with the following elements:
#' \itemize{
#'  \item graph - The igraph
#'  \item params - A list that contains the parameters of the graph generation
#'  and general graph information (node, edge numbers).
#' }
#'
#' @export
get_diffcor_graph <- S7::new_generic(
  name = "get_diffcor_graph",
  dispatch_args = "object",
  fun = function(object, min_cor = 0.2, fdr_threshold = 0.05, .verbose = TRUE) {
    S7::S7_dispatch()
  }
)

#' @export
#' @method get_diffcor_graph bulk_coexp
S7::method(get_diffcor_graph, bulk_coexp) <- function(
  object,
  min_cor = 0.2,
  fdr_threshold = 0.05,
  .verbose = TRUE
) {
  # Scope checks
  . <- delta_cor <- cor_a <- cor_b <- fdr <- weight <- NULL
  # Checks
  checkmate::assertClass(object, "bixverse::bulk_coexp")
  checkmate::qassert(min_cor, "R1[0,1]")
  checkmate::qassert(fdr_threshold, "R1[0,1]")
  checkmate::qassert(.verbose, "B1")
  # Function body
  cor_res <- S7::prop(object, "processed_data")[["differential_cor_res"]]
  graph_df <- cor_res$get_data_table(.verbose = .verbose) %>%
    .[, delta_cor := cor_a - cor_b] %>%
    .[, `:=`(
      cor_a = abs(cor_a),
      cor_b = abs(cor_b),
      fdr = rs_fdr_adjustment(p_val)
    )] %>%
    .[fdr <= fdr_threshold & (cor_a >= min_cor | cor_b >= min_cor)] %>%
    .[,
      weight := rs_range_norm(abs(delta_cor), max_val = 1, min_val = 0.05)
    ] %>%
    .[, c("feature_a", "feature_b")] %>%
    data.table::setnames(
      .,
      old = c("feature_a", "feature_b"),
      new = c("from", "to")
    )

  if (.verbose) {
    message(sprintf(
      "Generating differential correlation-based graph with %i edges.",
      nrow(graph_df)
    ))
  }

  graph <- igraph::graph_from_data_frame(graph_df, directed = FALSE)
  graph <- igraph::simplify(graph)

  graph_params <- list(
    min_cor = min_cor,
    fdr_threshold = fdr_threshold,
    no_nodes = length(igraph::V(graph)),
    no_edges = length(igraph::E(graph))
  )

  list(
    graph = graph,
    params = graph_params
  )
}

## coremo helpers --------------------------------------------------------------

#' Coremo: measures the quality of the clusters
#'
#' @description
#' Utility functions to give back the summary stats (median and mean R²) for
#' each of the identified clusters. Clusters over a size of 1000 genes will
#' be randomly sampled.
#'
#' @param modules Named vector. The names reflect the gene of the associated
#' module.
#' @param cor_mat Numeric matrix. The original correlation matrix.
#' @param random_seed Integer. Random seed to ensure consistency if sampling is
#' used.
#'
#' @return A data.table with the quality measures of the cluster.
coremo_cluster_quality <- function(modules, cor_mat, random_seed = 10101L) {
  # Checks
  checkmate::qassert(modules, c("S+", "I+"))
  checkmate::assertNamed(modules)
  checkmate::assertMatrix(cor_mat, mode = "numeric")
  checkmate::qassert(random_seed, "I1")
  # Function body
  cluster_list <- split(
    names(modules),
    modules
  )
  res <- rs_coremo_quality(
    cluster_genes = cluster_list,
    cor_mat = cor_mat,
    row_names = rownames(cor_mat),
    seed = random_seed
  ) %>%
    data.table::setDT()

  res
}

#' Coremo: cuts a hierarchical cluster based on k
#'
#' @description
#' This function uses a tree (output of [stats::hclust()]) and cuts it according
#' to the parameter k. If a `min_size` is specified, modules are merged by their
#' similarity of their eigen values in the distance matrix.
#'
#' @param tree hclust object. The hierarchical clustering of the correlation
#' matrix (or the distance thereof).
#' @param k Integer. Number of cuts on the tree.
#' @param dist_mat Numerical matrix. The distance matrix that was used to
#' compute the hierarchical clustering.
#' @param min_size Integer. Optional minimum size for the clusters.
#' @param cor_method String. Which correlation method to use for
#' optionally combining the small clusters. One of `c("pearson", "spearman")`.
#'
#' @return A vector with module membership.
#'
#' @importFrom magrittr `%>%`
coremo_tree_cut <- function(
  tree,
  k,
  dist_mat,
  min_size = NULL,
  cor_method = c("pearson", "spearman")
) {
  # Checks
  checkmate::assertClass(tree, "hclust")
  checkmate::qassert(k, "I1")
  checkmate::qassert(min_size, c("I1", "0"))
  checkmate::assertMatrix(dist_mat, mode = "numeric")
  checkmate::assertChoice(
    cor_method,
    c("pearson", "spearman")
  )
  # Function body
  clusters <- cutree(tree, k = k)
  # Early returns
  if (is.null(min_size)) {
    return(clusters)
  }
  cluster_size <- table(clusters)
  if (min(cluster_size) >= min_size) {
    return(clusters)
  }
  # Merge smaller clusters together
  to_keep <- names(cluster_size)[which(cluster_size >= min_size)]
  to_merge <- names(cluster_size)[which(cluster_size < min_size)]

  eg <- purrr::map(
    sort(unique(clusters)),
    \(cluster_i) {
      x <-
        prcomp(t(d[names(clusters)[clusters == cluster_i], ]), 1)$x[, "PC1"]
    }
  ) %>%
    do.call(rbind, .) %>%
    `rownames<-`(sort(unique(clusters)))

  spearman <- cor_method == "spearman"

  eg_cor <- rs_cor(x = t(eg), spearman = spearman)
  eg_cor <-
    abs(eg_cor[to_keep, toMerge, drop = FALSE])
  res <- clusters
  for (i in to_merge) {
    sel_clust <- to_keep[which.max(eg_cor[, as.character(i)])]
    res[which(res == as.numeric(i))] <- sel_clust
  }

  res
}

#' Coremo: Iterate over k for gene module detection.
#'
#' @description
#' This function uses a tree (output of [stats::hclust()]) and cuts it according
#' across all values from k_min to k_max and returns the quality at each
#' individual cut.
#'
#' @param tree hclust object. The hierarchical clustering of the correlation
#' matrix (or the distance thereof).
#' @param cor_mat Numerical matrix. Correlation matrix.
#' @param dist_mat Numerical matrix. Distance matrix.
#' @param k_min,k_max Integer. The minimum and maximum number of cuts.
#' @param min_size Integer. Optional minimum size of resulting modules.
#' @param cor_method String. Method for the correlation function. One of
#' `c("pearson", "spearman")`.
#' @param seed Integer. For reproducibility purposes.
#'
#' @return a data.table with stats (median size of the clusters, median weighted
#' R^2, and median R^2) on the varying levels of k.
#'
#' @importFrom magrittr `%>%`
#' @import data.table
tree_cut_iter <- function(
  tree,
  cor_mat,
  dist_mat,
  k_min = 1L,
  k_max = 100L,
  min_size = NULL,
  cor_method = c("spearman", "pearson"),
  seed = 42L
) {
  # Checks
  checkmate::assertClass(tree, "hclust")
  checkmate::assertMatrix(cor_mat, mode = "numeric")
  checkmate::assertMatrix(dist_mat, mode = "numeric")
  checkmate::qassert(k_min, "I1")
  checkmate::qassert(k_max, "I1")
  checkmate::qassert(min_size, c("I1", "0"))
  checkmate::assertChoice(cor_method, c("spearman", "pearson"))
  # Function body

  res <- purrr::map(
    k_min:k_max,
    \(k) {
      modules <-
        coremo_tree_cut(
          tree = tree,
          k = k,
          min_size = min_size,
          dist_mat = dist_mat,
          cor_method = cor_method
        ) %>%
        `names<-`(rownames(cor_mat))

      qc <- coremo_cluster_quality(
        modules = modules,
        cor_mat = cor_mat,
        random_seed = seed
      )

      res <- qc[, .(
        k = k,
        n = .N,
        "size_median" = median(size),
        "R2_weighted_median" = sum(r2med * size) / sum(size),
        "R2_median" = median(r2med)
      )]

      res
    }
  ) %>%
    data.table::rbindlist()

  res
}


#' Create distance object from a vector
#'
#' @param x Numerical vector. The upper-triangle values for which to generate
#' the `dist` object.
#' @param size Integer. Nrow (or ncol) of the symmetric matrix.
#'
#' @return Returns the distance object
create_dist_obj <- function(x, size) {
  # checks
  checkmate::qassert(x, "N+")
  checkmate::qassert(size, "I1")
  # body
  res <- x
  attr(res, "Size") <- size
  attr(res, "Diag") <- FALSE
  attr(res, "Upper") <- FALSE
  class(res) <- "dist"

  res
}

## getters ---------------------------------------------------------------------

#' @title Return the resolution results
#'
#' @description
#' Getter function to get the resolution results (if available).
#'
#' @param object The class, see [bixverse::bulk_coexp()].
#'
#' @return If resolution results were found, returns the data.table. Otherwise,
#' throws a warning and returns NULL.
#'
#' @export
get_resolution_res <- S7::new_generic(
  name = "get_resolution_res",
  dispatch_args = "object",
  fun = function(object) {
    S7::S7_dispatch()
  }
)

#' @export
#' @method get_resolution_res bulk_coexp
S7::method(get_resolution_res, bulk_coexp) <- function(object) {
  # Checks
  checkmate::assertClass(object, "bixverse::bulk_coexp")
  # Body
  resolution_results <- S7::prop(object, "outputs")[["resolution_results"]]
  if (is.null(resolution_results)) {
    warning(
      paste(
        "No resolution results found.",
        "Did you run cor_module_graph_check_res()? Returning NULL."
      )
    )
  }

  resolution_results
}

## plotting --------------------------------------------------------------------

# This one has a shared generic...

#' @export
#'
#' @import ggplot2
#'
#' @method plot_resolution_res bulk_coexp
S7::method(plot_resolution_res, bulk_coexp) <- function(
  object,
  print_head = TRUE,
  ...
) {
  # Checks
  checkmate::assertClass(object, "bixverse::bulk_coexp")
  checkmate::qassert(print_head, "B1")
  # Body
  plot_df <- S7::prop(object, "outputs")[["resolution_results"]]
  if (is.null(plot_df)) {
    warning(
      paste(
        "No resolution results found.",
        "Did you run cor_module_graph_check_res()? Returning NULL."
      )
    )
    return(NULL)
  }
  plot_df <- data.table::setorder(plot_df, -modularity)
  if (print_head) {
    print(head(plot_df))
  }
  p <- ggplot2::ggplot(
    data = plot_df,
    mapping = ggplot2::aes(x = resolution, y = modularity)
  ) +
    ggplot2::geom_point(
      mapping = ggplot2::aes(
        size = log10(good_clusters),
        fill = log10(avg_size)
      ),
      shape = 21,
      alpha = .7
    ) +
    ggplot2::xlab("Leiden cluster resolution") +
    ggplot2::ylab("Modularity") +
    ggplot2::theme_minimal() +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::scale_size_continuous(range = c(2, 8)) +
    ggplot2::labs(
      size = "Number of good clusters (log10)",
      fill = "Average cluster size (log10)"
    ) +
    ggplot2::ggtitle(
      "Resolution vs. modularity",
      subtitle = "With cluster number and size"
    )
  p
}


#' @title Plot the epsilon vs. power law goodness of fit result
#'
#' @description
#' Plots the epsilon results (if they can be found in the class). The x-axis
#' reflects the different epsilon parameters for the radial basis function,
#' and the y-axis the R2 value that the resulting networks follows a power law
#' distribution (i.e., scale free topology).
#'
#' @param object The class, see [bixverse::bulk_coexp()].
#'
#' @return If epsilon results were found, returns the ggplot. Otherwise, throws
#' a warning and returns NULL.
#'
#' @export
plot_epsilon_res <- S7::new_generic(
  name = "plot_epsilon_res",
  dispatch_args = "object",
  fun = function(object) {
    S7::S7_dispatch()
  }
)

#' @export
#'
#' @import ggplot2
#'
#' @method plot_epsilon_res bulk_coexp
S7::method(plot_epsilon_res, bulk_coexp) <- function(object) {
  # Checks
  checkmate::assertClass(object, "bixverse::bulk_coexp")
  # Body
  plot_df <- S7::prop(object, "outputs")[["epsilon_data"]]
  if (is.null(plot_df)) {
    warning(
      paste(
        "No resolution results found.",
        "Did you run cor_module_check_epsilon()? Returning NULL."
      )
    )
    return(NULL)
  }

  p <- ggplot2::ggplot(data = plot_df, ggplot2::aes(x = epsilon, y = r2_vals)) +
    ggplot2::geom_point(size = 3, shape = 21) +
    ggplot2::geom_line() +
    ggplot2::theme_minimal() +
    ggplot2::ylim(0, 1) +
    ggplot2::xlab("Epsilon") +
    ggplot2::ylab("Goodness of fit (R2)") +
    ggplot2::ggtitle(
      "Epsilon vs. scale free topology",
      subtitle = "Goodness of fit for log(connectivity) ~ log(p(connectivity)"
    )

  p
}


#' @title Plot the k cuts vs median R2
#'
#' @description
#' Plots the optimal k vs. median of median R2 graph to identify the optimal
#' number of cuts.
#'
#' @param object The class, see [bixverse::bulk_coexp()].
#'
#' @return If optimal cuts results were found, returns the ggplot. Otherwise,
#' throws a warning and returns NULL.
#'
#' @export
plot_optimal_cuts <- S7::new_generic(
  name = "plot_optimal_cuts",
  dispatch_args = "object",
  fun = function(object) {
    S7::S7_dispatch()
  }
)

#' @export
#'
#' @import ggplot2
#'
#' @method plot_optimal_cuts bulk_coexp
S7::method(plot_optimal_cuts, bulk_coexp) <- function(object) {
  # Checks
  checkmate::assertClass(object, "bixverse::bulk_coexp")

  plot_df <- data.table::copy(S7::prop(object = object, name = "outputs")[[
    "optimal_cuts"
  ]]) %>%
    data.table::setorder(gradient_change)

  if (is.null(plot_df)) {
    warning(paste(
      "No optimal_cuts data.table found.",
      "Did you run cor_module_coremo_clustering()?",
      "Returning NULL."
    ))
  }

  optimal_cuts <- S7::prop(object = object, name = "params")[["coremo"]][[
    "inflection_idx"
  ]]

  p <- ggplot2::ggplot(
    data = plot_df,
    mapping = ggplot2::aes(x = k, y = R2_weighted_median)
  ) +
    ggplot2::geom_point(
      mapping = ggplot2::aes(fill = gradient_change),
      shape = 21,
      alpha = 0.7,
      size = 3
    ) +
    ggplot2::scale_fill_viridis_c() +
    ggplot2::theme_minimal() +
    ggplot2::xlab("k cuts") +
    ggplot2::ylab("Median of median weighted R2") +
    ggplot2::labs(fill = "Gradient change") +
    ggplot2::geom_vline(
      xintercept = optimal_cuts,
      color = "darkgrey",
      linetype = "dashed"
    ) +
    ggplot2::ggtitle("k cuts vs. change in R2")

  p
}
