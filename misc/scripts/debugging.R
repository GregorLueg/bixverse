object = cor_test
graph_params = params_cor_graph(epsilon = 2.5)
resolution_params = params_graph_resolution()
random_seed = 123L
min_genes = 10L
parallel = TRUE
max_workers = NULL
.verbose = TRUE

detection_method <- S7::prop(object, "params")[["detection_method"]]


c(graph, graph_params) %<-%
  with(
    graph_params,
    switch(
      detection_method,
      "correlation-based" = get_cor_graph(
        object = object,
        epsilon = epsilon,
        .verbose = .verbose
      ),
      "differential correlation-based" = get_diffcor_graph(
        object = object,
        min_cor = min_cor,
        fdr_threshold = fdr_threshold,
        .verbose = .verbose
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
    message(sprintf(
      "Using parallel computation over %i cores via mirai.",
      max_workers
    ))
  }

  mirai::daemons(max_workers)

  community_df_res <- mirai::mirai_map(
    resolutions,
    \(res, seed, graph) {
      set.seed(seed)
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
    .args = list(seed = random_seed, graph = graph)
  )[]

  mirai::daemons(0)
} else {
  if (.verbose) {
    message("Using sequential computation.")
  }
  community_df_res <- purrr::map(
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
    }
  )
}
