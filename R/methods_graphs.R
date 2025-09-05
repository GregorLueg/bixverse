# network_diffusions -----------------------------------------------------------

## diffusion methods -----------------------------------------------------------

#' Diffuse seed genes over a network
#'
#' @description
#' This function takes a diffusion vector and leverages personalised page-rank
#' diffusion to identify influential nodes. These can be used subsequently for
#' community detection or check AUROC values given a set of genes.
#'
#' @param object `network_diffusions` object. The underlying class
#' [bixverse::network_diffusions()].
#' @param diffusion_vector Named nuermic. A named vector with values to use for
#' the reset parameter in the personalised page-rank diffusion. Names should
#' represent node names of the graph.
#' @param summarisation String. If there are duplicated names in the
#' `diffusion_vector` how to summarise the scores.
#'
#' @return The class with added diffusion score based on a single set of seed
#' genes. Additionally, the seed genes are stored in the class.
#'
#' @export
diffuse_seed_nodes <- S7::new_generic(
  name = "diffuse_seed_nodes",
  dispatch_args = "object",
  fun = function(
    object,
    diffusion_vector,
    summarisation = c("max", "mean", "harmonic_sum")
  ) {
    S7::S7_dispatch()
  }
)

#' @export
#'
#' @importFrom magrittr `%>%`
#'
#' @method diffuse_seed_nodes network_diffusions
S7::method(diffuse_seed_nodes, network_diffusions) <-
  function(
    object,
    diffusion_vector,
    summarisation = c("max", "mean", "harmonic_sum")
  ) {
    # Checks
    checkmate::assertClass(object, "bixverse::network_diffusions")
    checkmate::assertNumeric(diffusion_vector)
    checkmate::assertNamed(diffusion_vector, .var.name = "diffusion_vector")
    checkmate::assertChoice(summarisation, c("max", "mean", "harmonic_sum"))

    # Body
    ## Create the diffusion vector
    diffusion_vector <- summarise_scores(
      diffusion_vector,
      summarisation = summarisation
    )
    nodes_names <- igraph::V(S7::prop(object, "graph"))$name
    seed_nodes <- intersect(names(diffusion_vector), nodes_names)
    diff_vec <- rep(0, length(nodes_names)) %>% `names<-`(nodes_names)
    for (node in seed_nodes) {
      diff_vec[node] <- diffusion_vector[node]
    }

    if (sum(diff_vec) == 0) {
      stop(
        paste(
          "No scores found to diffuse over the network.",
          "Please check the names and/or values."
        )
      )
    }

    ## Create the page-rank diffusion
    page_rank_score <- igraph::page_rank(
      S7::prop(object, "graph"),
      personalized = diff_vec
    )

    diffusion_params <- list(
      "seed_nodes" = seed_nodes,
      "diffusion_type" = "single",
      "diffusion_vector" = diffusion_vector
    )

    ## Assign and return
    S7::prop(object, "diffusion_res") <- page_rank_score$vector
    S7::prop(object, "params")[["diffusion_params"]] <- diffusion_params

    return(object)
  }


#' Diffuse seed genes in a tied manner over a network
#'
#' @description
#' This function takes two sets of diffusion vector and leverages tied diffusion
#' to identify an intersection of influential nodes. If the network is
#' undirected, the method will run two personalised page rank diffusions based
#' on the diffusion vectors and generate the score aggregation
#'
#' @param object `network_diffusions` object. The underlying class
#' [bixverse::network_diffusions()].
#' @param diffusion_vector_1 Named numeric. The first named vector with values
#' to use for the reset parameter in the personalised page-rank diffusion. Names
#' should represent node names of the graph.
#' @param diffusion_vector_2 Named numeric. The second named vector with values
#' to use for the reset parameter in the personalised page-rank diffusion. Names
#' should represent node names of the graph.
#' @param summarisation String. If there are duplicated names in the
#' `diffusion_vector` how to summarise these.
#' @param score_aggregation String. How to summarise the tied scores.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @return The class with added diffusion score based on a two sets of seed
#' genes. Additionally, the seed genes are stored in the class.
#'
#' @export
tied_diffusion <- S7::new_generic(
  name = "tied_diffusion",
  dispatch_args = "object",
  fun = function(
    object,
    diffusion_vector_1,
    diffusion_vector_2,
    summarisation = c("max", "mean", "harmonic_sum"),
    score_aggregation = c("min", "max", "mean"),
    .verbose = FALSE
  ) {
    S7::S7_dispatch()
  }
)

#' @export
#'
#' @importFrom magrittr `%>%`
#'
#' @method tied_diffusion network_diffusions
S7::method(tied_diffusion, network_diffusions) <-
  function(
    object,
    diffusion_vector_1,
    diffusion_vector_2,
    summarisation = c("max", "mean", "harmonic_sum"),
    score_aggregation = c("min", "max", "mean"),
    .verbose = FALSE
  ) {
    # Checks
    checkmate::assertClass(object, "bixverse::network_diffusions")
    checkmate::assertNumeric(diffusion_vector_1)
    checkmate::assertNamed(diffusion_vector_1, .var.name = "diffusion_vector_1")
    checkmate::assertNumeric(diffusion_vector_2)
    checkmate::assertNamed(diffusion_vector_2, .var.name = "diffusion_vector_2")
    checkmate::assertChoice(summarisation, c("max", "mean", "harmonic_sum"))
    checkmate::assertChoice(score_aggregation, c("min", "max", "mean"))
    checkmate::qassert(.verbose, "B1")

    # Body
    ## Create the diffusion vectors
    diffusion_vector_1 <- summarise_scores(
      diffusion_vector_1,
      summarisation = summarisation
    )
    diffusion_vector_2 <- summarise_scores(
      diffusion_vector_2,
      summarisation = summarisation
    )
    nodes_names <- igraph::V(S7::prop(object, "graph"))$name
    seed_nodes_1 <- intersect(names(diffusion_vector_1), nodes_names)
    seed_nodes_2 <- intersect(names(diffusion_vector_2), nodes_names)
    diff_vec_1 <- diff_vec_2 <- rep(0, length(nodes_names)) %>%
      `names<-`(nodes_names)
    for (node in seed_nodes_1) {
      diff_vec_1[node] <- diffusion_vector_1[node]
    }
    for (node in seed_nodes_2) {
      diff_vec_2[node] <- diffusion_vector_2[node]
    }
    if ((sum(diff_vec_1) == 0) || (sum(diff_vec_1) == 0)) {
      stop(
        paste(
          "No scores found on first and/or second of the diffusion vectors.",
          "Please check the names and/or values."
        )
      )
    }

    ## First diffusion
    score_1 <- igraph::page_rank(
      S7::prop(object, "graph"),
      personalized = diff_vec_1
    )$vector

    ## Second diffusion
    directed <- S7::prop(object, "params")[["directed_graph"]]
    score_2 <- if (directed) {
      if (.verbose) {
        message(
          paste(
            "Directed graph found.",
            "Function will use transpose of adjacency for second diffusion."
          )
        )
      }
      adj <- igraph::as_adjacency_matrix(S7::prop(object, "graph"))
      adj_t <- Matrix::t(adj)
      igraph_obj_t <- igraph::graph_from_adjacency_matrix(adj_t)
      igraph::page_rank(igraph_obj_t, personalized = diff_vec_2)$vector
    } else {
      if (.verbose) {
        message(
          paste(
            "Undirected graph found.",
            "Using graph as is for second diffusion."
          )
        )
      }
      igraph::page_rank(
        S7::prop(object, "graph"),
        personalized = diff_vec_2
      )$vector
    }

    ## Summarise the scores
    final_tiedie_diffusion <- switch(
      score_aggregation,
      "min" = pmin(score_1, score_2),
      "max" = pmax(score_1, score_2),
      rowMeans(cbind(score_1, score_2))
    )

    diffusion_params <- list(
      "seed_nodes_1" = seed_nodes_1,
      "seed_nodes_2" = seed_nodes_2,
      "diffusion_type" = "tied",
      "diffusion_vector_1" = diffusion_vector_1,
      "diffusion_vector_2" = diffusion_vector_2,
      "score_aggregation" = score_aggregation
    )

    ## Assign and return
    S7::prop(object, "diffusion_res") <- final_tiedie_diffusion
    S7::prop(object, "params")[[
      "diffusion_params"
    ]] <- diffusion_params

    return(object)
  }


#' Generate permuation scores for the diffusion
#'
#' @description
#' This function generate node-degree adjusted permutations of a given diffusion
#' score and adds Z-scores to the object. The function will automatically
#' determine if the original diffusion was a single or tied diffusion and
#' construct permutations accordingly.
#'
#' @param object `network_diffusions` object. The underlying class
#' [bixverse::network_diffusions()].
#' @param perm_iters Integer. Number of permutations to test for. Defaults to
#' `1000L`.
#' @param random_seed Integer. Random seed for determinism.
#' @param .verbose Boolean. Controls verbosity.
#'
#' @return The class with added diffusion score based on a single set of seed
#' genes. Additionally, the seed genes are stored in the class.
#'
#' @export
permute_seed_nodes <- S7::new_generic(
  name = "permute_seed_nodes",
  dispatch_args = "object",
  fun = function(
    object,
    perm_iters = 1000L,
    random_seed = 10101L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)

#' @export
#'
#' @importFrom magrittr `%>%`
#'
#' @method permute_seed_nodes network_diffusions
S7::method(permute_seed_nodes, network_diffusions) <- function(
  object,
  perm_iters = 1000L,
  random_seed = 10101L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertClass(object, "bixverse::network_diffusions")
  checkmate::qassert(perm_iters, "I1")
  checkmate::qassert(random_seed, "I1")
  checkmate::qassert(.verbose, "B1")

  # function body
  graph <- S7::prop(object, "graph")
  diffusion_params <- S7::prop(object, "params")[["diffusion_params"]]
  diffusion_results <- S7::prop(object, "diffusion_res")
  nodes_names <- igraph::V(graph)$name

  # prepare data for rust
  edge_list <- igraph::as_edgelist(graph, names = TRUE)
  graph_names <- igraph::V(graph)$name

  weights <- if (igraph::is_weighted(graph)) {
    igraph::edge.attributes(graph)$weight
  } else {
    NULL
  }

  if (diffusion_params$diffusion_type == "single") {
    if (.verbose) {
      message(sprintf(
        "Permutations for single diffusion with %i iters being generated.",
        perm_iters
      ))
    }

    # generate randomised diffusion vecs
    diffusion_vector <- diffusion_params[["diffusion_vector"]]

    randomised_diffusions <- generate_perm_diffusion_vecs(
      graph = graph,
      diffusion_vec = diffusion_vector,
      iters = perm_iters
    )

    # use rust for fast calculations
    page_rank_perm_res <- rs_page_rank_parallel(
      node_names = graph_names,
      from = edge_list[, 1],
      to = edge_list[, 2],
      weights = weights,
      diffusion_scores = randomised_diffusions,
      undirected = !igraph::is_directed(graph)
    )
  } else {
    if (.verbose) {
      message(sprintf(
        "Permutations for tied diffusion with %i iters being generated.",
        perm_iters
      ))
    }
    diffusion_vector_1 <- diffusion_params[["diffusion_vector_1"]]
    diffusion_vector_2 <- diffusion_params[["diffusion_vector_2"]]

    permutations_1 <- generate_perm_diffusion_vecs(
      graph = graph,
      diffusion_vec = diffusion_vector_1,
      iters = perm_iters
    )

    permutations_2 <- generate_perm_diffusion_vecs(
      graph = graph,
      diffusion_vec = diffusion_vector_2,
      iters = perm_iters
    )

    page_rank_perm_res <- rs_tied_diffusion_parallel(
      node_names = graph_names,
      from = edge_list[, 1],
      to = edge_list[, 2],
      weights = weights,
      diffusion_scores_1 = permutations_1,
      diffusion_scores_2 = permutations_2,
      summarisation_fun = diffusion_params$score_aggregation,
      undirected = !igraph::is_directed(graph)
    )
  }

  diffuion_means <- colMeans(page_rank_perm_res)
  diffusion_sds <- matrixStats::colSds(page_rank_perm_res)

  z_scores <- (diffusion_results - diffuion_means) /
    (diffusion_sds + 10^-32)

  diffusion_perm_params <- list(
    "perm_iters" = perm_iters,
    "random_seed" = random_seed,
    "perm_mean" = diffuion_means,
    "perm_sds" = diffusion_sds
  )

  S7::prop(object, "diffusion_perm") <- z_scores
  S7::prop(object, "params")[["diffusion_perm_params"]] <- diffusion_perm_params

  return(object)
}


## community detection ---------------------------------------------------------

#' Identify privileged communities based on a given diffusion vector
#'
#' @description Detects privileged communities after a diffusion based on seed
#' nodes.
#'
#' @param object `network_diffusions` object. The underlying class
#' [bixverse::network_diffusions()].
#' @param community_params List. Parameters for the community detection within
#' the reduced network, see [bixverse::params_community_detection()]. A list
#' with the following items:
#' \itemize{
#'  \item max_nodes - Integer. Number of maximum nodes per community. Larger
#'  communities will be recursively subclustered.
#'  \item min_nodes - Integer. Minimum number of nodes per community.
#'  \item min_seed_nodes - Integer. Minimum number of seed genes that have to
#'  be found in a given community.
#'  \item initial_res - Float. Initial resolution parameter for the Leiden
#'  clustering.
#'  \item threshold_type - String. One of `c("prop_based", "pval_based")`.
#'  You can chose to include a certain proportion of the network (like in the
#'  original paper) with the highest diffusion scores, or use p-values based
#'  on permutations. Defaults to `"prop_based"`.
#'  \item network_threshold - Float. The proportion of the network to
#'  include. Used if `threshold_type = "prop_based"`.
#'  \item pval_threshold - Float. The maximum p-value for nodes to be included.
#'  Used if `threshold_type = "pval_based"`.
#' }
#' @param seed Random seed.
#' @param .verbose Controls the verbosity of the function.
#' @param .max_iters Controls how many iterations shall be tried for the
#' sub-clustering. To note, in each iteration of the sub-clustering, the
#' resolution parameter is increased by 0.05, to identify more granular
#' communities within the sub communities.
#'
#' @return The class with added diffusion community detection results (if any
#' could be identified with the provided parameters).
#'
#' @export
community_detection <- S7::new_generic(
  name = "community_detection",
  dispatch_args = "object",
  fun = function(
    object,
    community_params = params_community_detection(),
    seed = 42L,
    .verbose = FALSE,
    .max_iters = 100L
  ) {
    S7::S7_dispatch()
  }
)


#' @export
#'
#' @import data.table
#' @importFrom magrittr `%>%`
#'
#' @method community_detection network_diffusions
S7::method(community_detection, network_diffusions) <- function(
  object,
  community_params = params_community_detection(),
  seed = 42L,
  .verbose = FALSE,
  .max_iters = 100L
) {
  # Bindings
  . <- N <- cluster_id <- node_id <- cluster_size <- seed_nodes_no <-
    seed_nodes_no <- seed_nodes_1 <- seed_nodes_2 <- seed_node <- `:=` <- NULL
  # Checks
  checkmate::assertClass(object, "bixverse::network_diffusions")
  assertCommunityParams(community_params)
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, "B1")
  checkmate::qassert(.max_iters, "I1")

  # reduce the graph
  diffusion_score <- S7::prop(object, "diffusion_res")
  # early return
  if (length(diffusion_score) == 0) {
    warning(
      paste(
        "The diffusion score has length 0.",
        "Likely you did not run the diffusion methods. Returning class as is."
      )
    )
    return(object)
  }

  # subset the graph
  if (community_params$threshold_type == "prop_based") {
    nodes_to_include <- diffusion_score %>%
      sort(decreasing = TRUE) %>%
      .[1:ceiling(community_params$network_threshold * length(diffusion_score))]
  } else {
    z_score <- S7::prop(object, "diffusion_perm")
    if (length(z_score) == 0) {
      warning("No z-scores found. Will calculate these now.")
      object <- permute_seed_nodes(object)
      z_score <- S7::prop(object, "diffusion_perm")
    }
    p_vals <- pnorm(q = abs(z_score), lower.tail = FALSE)
    to_include <- p_vals <= community_params$pval_threshold
    nodes_to_include <- diffusion_score[to_include]
  }

  red_graph <- igraph::subgraph(
    S7::prop(object, "graph"),
    names(nodes_to_include)
  )

  ## First clustering
  set.seed(seed)

  if (S7::prop(object, "params")$directed_graph) {
    red_graph <- igraph::as_undirected(red_graph, mode = "each")
  }

  final_clusters <- with(community_params, {
    first_clusters <- igraph::cluster_leiden(
      red_graph,
      resolution = initial_res,
      n_iterations = 5,
      objective_function = "modularity"
    )

    clusters_df <- data.table::data.table(
      node_id = first_clusters$names,
      cluster_id = first_clusters$membership
    )

    node_frequency <- clusters_df[, .N, .(cluster_id)]

    ## Subclustering
    clusters_with_too_many_nodes <- node_frequency[N > max_nodes, cluster_id]
    final_clusters <- clusters_df[!cluster_id %in% clusters_with_too_many_nodes]

    for (i in seq_along(clusters_with_too_many_nodes)) {
      cluster_i <- clusters_with_too_many_nodes[i]
      nodes_in_cluster <- clusters_df[cluster_id == cluster_i, node_id]
      finalised_clusters <- data.table()
      # Loop through, until all clusters are below the minimum genes or max
      # iterations is hit
      l <- 1
      while (length(nodes_in_cluster) != 0) {
        set.seed(seed + l)
        if (.verbose) {
          message("Cluster ", i, " gets subclustered. Iter: ", l)
        }
        red_graph_l <- igraph::subgraph(
          red_graph,
          data.table::chmatch(nodes_in_cluster, igraph::V(red_graph)$name)
        )

        clusters_red <- igraph::cluster_leiden(
          red_graph_l,
          resolution = initial_res + l * 0.05,
          n_iterations = 5,
          objective_function = "modularity"
        )

        subclusters <- data.table(
          node_id = clusters_red$names,
          cluster_id = clusters_red$membership
        )
        subclusters_frequency <- subclusters[, .N, .(cluster_id)]
        clusters_small_enough <- subclusters_frequency[
          N <= max_nodes,
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

    return(final_clusters)
  })

  ## Add the seed node information based on diffusion type
  diffusion_params <- S7::prop(object, "params")[[
    "diffusion_params"
  ]]

  diffusion_type <- diffusion_params$diffusion_type

  final_node_frequency <- with(community_params, {
    if (diffusion_type == "single") {
      seed_nodes <- diffusion_params$seed_nodes

      final_clusters[,
        .(
          cluster_size = length(node_id),
          seed_nodes_no = sum(node_id %in% seed_nodes)
        ),
        .(cluster_id)
      ]
    } else {
      seed_nodes_set_1 <- diffusion_params$seed_nodes_1
      seed_nodes_set_2 <- diffusion_params$seed_nodes_2

      final_clusters[,
        .(
          cluster_size = length(node_id),
          seed_nodes_1 = sum(node_id %in% seed_nodes_set_1),
          seed_nodes_2 = sum(node_id %in% seed_nodes_set_2)
        ),
        .(cluster_id)
      ]
    }
  })

  ## Finalise the clusters
  clusters_to_take <- with(community_params, {
    if (diffusion_type == "single") {
      final_node_frequency[
        cluster_size >= min_nodes &
          seed_nodes_no >= min_seed_nodes,
        cluster_id
      ]
    } else {
      final_node_frequency[
        cluster_size >= min_nodes &
          seed_nodes_1 >= min_seed_nodes &
          seed_nodes_2 >= min_seed_nodes,
        cluster_id
      ]
    }
  })

  # Early return
  if (length(clusters_to_take) == 0) {
    warning(
      paste(
        "No communities found with the given parameters.",
        "Returning class as is."
      )
    )
    return(object)
  }

  finalised_clusters_clean <- final_clusters[cluster_id %in% clusters_to_take]

  ks_vals <- vector(mode = "numeric", length = length(clusters_to_take))

  for (i in seq_along(clusters_to_take)) {
    cluster <- clusters_to_take[i]
    cluster_nodes <- finalised_clusters_clean[cluster_id == cluster, node_id]
    ks <- suppressWarnings(ks.test(
      diffusion_score[cluster_nodes],
      diffusion_score[which(!names(diffusion_score) %in% cluster_nodes)],
      alternative = "less"
    ))
    ks_vals[i] <- ks$p.value
  }

  ks_val_df <- data.table::data.table(
    cluster_id = clusters_to_take,
    ks_pval = ks_vals
  )

  final_result <- purrr::reduce(
    list(finalised_clusters_clean, ks_val_df, final_node_frequency),
    merge,
    by = "cluster_id"
  ) %>%
    .[, `:=`(
      diffusion_score = diffusion_score[node_id],
      cluster_id = as.character(cluster_id)
    )]

  if (diffusion_type == "single") {
    seed_nodes <- diffusion_params$seed_nodes
    final_result[, seed_node := node_id %in% seed_nodes]
  } else {
    seed_nodes_set_1 <- diffusion_params$seed_nodes_1
    seed_nodes_set_2 <- diffusion_params$seed_nodes_2
    final_result[, `:=`(
      seed_node_a = node_id %in% seed_nodes_set_1,
      seed_node_b = node_id %in% seed_nodes_set_2
    )]
  }

  cluster_name_prettifier <- setNames(
    paste(
      "cluster",
      seq_along(unique(
        final_result$cluster_id
      )),
      sep = "_"
    ),
    unique(final_result$cluster_id)
  )

  final_result[, cluster_id := cluster_name_prettifier[cluster_id]]

  ## Assign and return
  S7::prop(object, "final_results") <- final_result
  S7::prop(object, "params")[["community_params"]] <- with(
    community_params,
    list(
      network_threshold = network_threshold,
      pval_threshold = pval_threshold,
      threshold_type = threshold_type,
      max_nodes = max_nodes,
      min_nodes = min_seed_nodes,
      min_seed_nodes = min_seed_nodes
    )
  )

  return(object)
}

## utils -----------------------------------------------------------------------

#' Calculate the AUROC for a diffusion score
#'
#' @description
#' This functions can take a given `network_diffusions` object and calculates an
#' AUC and generates a Z-score based on random permutation of `random_aucs` for
#' test for statistical significance if desired.
#'
#' @param object `network_diffusions` object. The underlying class
#' [bixverse::network_diffusions()].
#' @param hit_nodes String vector. Which nodes in the graph are considered a
#' 'hit'.
#' @param auc_iters Integer. How many iterations to run to approximate the
#' AUROC.
#' @param random_aucs Integer. How many random AUROCs to calculate to estimate
#' the Z-score. Only of relevance if permutation test is set to `TRUE`.
#' @param permutation_test Boolean. Shall a permutation based Z-score be
#' calculated.
#' @param seed Integer. Random seed.
#'
#' @return List with AUC and Z-score as the two named elements if permutations
#' test set to TRUE; otherwise just the AUC.
#'
#' @export
calculate_diffusion_auc <- S7::new_generic(
  name = "calculate_diffusion_auc",
  dispatch_args = "object",
  fun = function(
    object,
    hit_nodes,
    auc_iters = 10000L,
    random_aucs = 1000L,
    permutation_test = FALSE,
    seed = 42L
  ) {
    S7::S7_dispatch()
  }
)


#' @export
#'
#' @importFrom magrittr `%>%`
#'
#' @method tied_diffusion network_diffusions
S7::method(calculate_diffusion_auc, network_diffusions) <-
  function(
    object,
    hit_nodes,
    auc_iters = 10000L,
    random_aucs = 1000L,
    permutation_test = FALSE,
    seed = 42L
  ) {
    # Checks
    checkmate::assertClass(object, "bixverse::network_diffusions")
    checkmate::qassert(hit_nodes, "S+")
    checkmate::qassert(auc_iters, "I1")
    checkmate::qassert(seed, "I1")
    checkmate::qassert(permutation_test, "B1")
    if (permutation_test) {
      checkmate::qassert(random_aucs, "I1")
    }
    # Body
    diffusion_score <- S7::prop(object, "diffusion_res")
    if (length(diffusion_score) == 0) {
      warning(
        paste(
          "The diffusion score has length 0.",
          "Likely you did not run the diffusion methods. Returning NULL."
        )
      )
      return(NULL)
    }
    pos_scores <- diffusion_score[hit_nodes]
    neg_scores <- diffusion_score[which(!names(diffusion_score) %in% hit_nodes)]
    auc <- rs_fast_auc(
      pos_scores = pos_scores,
      neg_scores = neg_scores,
      iters = auc_iters,
      seed = seed
    )
    if (permutation_test) {
      random_aucs <- rs_create_random_aucs(
        score_vec = diffusion_score,
        size_pos = length(pos_scores),
        random_iters = random_aucs,
        auc_iters = auc_iters,
        seed = seed
      )

      z <- (auc - mean(random_aucs)) / sd(random_aucs)

      to_ret <- list(auc = auc, z = z)
    } else {
      to_ret <- auc
    }

    return(to_ret)
  }

### helpers --------------------------------------------------------------------

##### random page-rank perms ---------------------------------------------------

#' Generates random permutation vectors
#'
#' @param graph igraph. The graph for which to generate the random diffusion
#' vectors
#' @param diffusion_vec Named numeric. The initial diffusion vector.
#' @param bins Integer. Number of bins to use for node degree aware sampling.
#' @param iters Integer. Number of random permutations to generate.
#' @param random_seed Integer. Random seed.
#'
#' @return List with the permutations.
#'
#' @importFrom magrittr `%$%`
generate_perm_diffusion_vecs <- function(
  graph,
  diffusion_vec,
  bins = 25L,
  iters = 1000L,
  random_seed = 10101L
) {
  # checks
  checkmate::assertClass(graph, "igraph")
  checkmate::assertNumeric(diffusion_vec, names = "named")
  checkmate::qassert(bins, "I1")
  checkmate::qassert(iters, "I1")
  checkmate::qassert(random_seed, "I1")

  nodes_names <- igraph::V(graph)$name
  node_degree_distribution <- log(igraph::degree(graph))
  node_degree_discrete <- cut(node_degree_distribution, bins) %>%
    `names<-`(names(node_degree_distribution))

  degree_groups <- split(names(node_degree_discrete), node_degree_discrete)
  diffusion_names <- intersect(names(diffusion_vec), nodes_names)
  node_degrees <- node_degree_discrete[diffusion_names]

  randomised_diffusions <- purrr::map(1:iters, function(i) {
    set.seed(random_seed + i)
    random_set_i <- vapply(
      node_degrees,
      function(degree) {
        sample(degree_groups[[as.character(degree)]], 1)
      },
      character(1)
    )

    diff_vec_i <- rep(0, length(nodes_names)) %>% `names<-`(nodes_names)

    diff_vec_i[random_set_i] <- diffusion_vec[names(random_set_i)]

    diff_vec_i / sum(diff_vec_i)
  })

  return(randomised_diffusions)
}

#### utils ---------------------------------------------------------------------

#' Summarise gene scores if they are duplicates.
#'
#' @param x Named numeric.
#' @param summarisation String. Which summary function to use.
#'
#' @return Named numeric.
#'
#' @importFrom magrittr `%$%`
summarise_scores <- function(
  x,
  summarisation = c("max", "mean", "harmonic_sum")
) {
  # devtools::check() stuff
  value <- . <- node_name <- setNames <- NULL
  # Checks
  checkmate::assertNumeric(x)
  checkmate::assertNamed(x, .var.name = "x")
  checkmate::assertChoice(summarisation, c("max", "mean", "harmonic_sum"))
  # Body
  dt <- data.table::data.table(node_name = names(x), value = x)
  summary_fun <- switch(
    summarisation,
    "mean" = rlang::expr(mean(value)),
    "max" = rlang::expr(max(value)),
    rlang::expr(bixverse::ot_harmonic_score(value)) # Default case
  )
  res <-
    rlang::eval_tidy(rlang::quo(dt[,
      .(value = !!summary_fun),
      .(node_name)
    ])) %$%
    setNames(value, node_name)
  res
}

# rbh_graph --------------------------------------------------------------------

## graph generation ------------------------------------------------------------

#' Generate an RBH graph.
#'
#' @description
#' This function will generate an RBH graph based on set similarity between
#' gene modules. You have the option to use an overlap coefficient instead of
#' Jaccard similarity and to specify a minimum similarity.
#'
#' @param object The underlying class, see [bixverse::rbh_graph()].
#' @param minimum_similarity The minimum similarity to create an edge.
#' @param overlap_coefficient Boolean. Shall the overlap coefficient be used
#' instead of Jaccard similarity. Only relevant if the underlying class is set
#' to set similarity.
#' @param spearman Boolean. Shall Spearman correlation be used. Only relevant
#' if the underlying class is set to correlation-based similarity.
#'
#' @return The class with added properties.
#'
#' @export
generate_rbh_graph <- S7::new_generic(
  name = "generate_rbh_graph",
  dispatch_args = "object",
  fun = function(
    object,
    minimum_similarity,
    overlap_coefficient = FALSE,
    spearman = FALSE
  ) {
    S7::S7_dispatch()
  }
)


#' @export
#'
#' @importFrom magrittr `%>%`
#'
#' @method generate_rbh_graph rbh_graph
S7::method(generate_rbh_graph, rbh_graph) <-
  function(
    object,
    minimum_similarity,
    overlap_coefficient = FALSE,
    spearman = FALSE
  ) {
    # scope checks
    origin_modules <- . <- similiarity <- origin <- target <- `:=` <-
      target_modules <- NULL
    # checks
    checkmate::assertClass(object, "bixverse::rbh_graph")
    checkmate::qassert(minimum_similarity, "R[0, 1]")
    checkmate::qassert(overlap_coefficient, "B1")
    checkmate::qassert(spearman, "B1")

    # function body
    rbh_type <- S7::prop(object, "params")[["rbh_type"]]

    rbh_results <- if (rbh_type == "set") {
      list_of_list <- S7::prop(object, "module_data")

      rs_rbh_sets(
        module_list = list_of_list,
        overlap_coefficient = overlap_coefficient,
        min_similarity = minimum_similarity
      )
    } else {
      list_of_matrices <- S7::prop(object, "module_data")

      rs_rbh_cor(
        module_matrices = list_of_matrices,
        spearman = spearman,
        min_similarity = minimum_similarity
      )
    }

    rbh_results$origin_modules[rbh_results$origin_modules == "NA"] <- NA
    rbh_results$target_modules[rbh_results$target_modules == "NA"] <- NA

    origin_vector <- unlist(purrr::map2(
      rbh_results$origin,
      rbh_results$comparisons,
      ~ {
        rep(.x, each = .y)
      }
    ))

    target_vector <- unlist(purrr::map2(
      rbh_results$target,
      rbh_results$comparisons,
      ~ {
        rep(.x, each = .y)
      }
    ))

    rbh_results_dt <- data.table::data.table(
      origin = origin_vector,
      target = target_vector,
      origin_modules = rbh_results$origin_modules,
      target_modules = rbh_results$target_modules,
      similiarity = rbh_results$similarity
    ) %>%
      .[
        !is.na(origin_modules) &
          similiarity >= minimum_similarity
      ] %>%
      .[, `:=`(
        combined_origin = paste(origin, origin_modules, sep = "_"),
        combined_target = paste(target, target_modules, sep = "_")
      )]

    edge_dt <- rbh_results_dt[, c(
      "combined_origin",
      "combined_target",
      "similiarity"
    )] %>%
      data.table::setnames(
        .,
        old = c("combined_origin", "combined_target", "similiarity"),
        new = c("from", "to", "weight")
      )

    rbh_igraph <- igraph::graph_from_data_frame(edge_dt, directed = FALSE)

    ## Assign and return
    S7::prop(object, "rbh_edge_df") <- rbh_results_dt
    S7::prop(object, "rbh_graph") <- rbh_igraph
    S7::prop(object, "params")[["rbh_graph_gen"]] <- list(
      minimum_similarity = minimum_similarity,
      overlap_coefficient = overlap_coefficient,
      spearman = spearman
    )

    return(object)
  }


#' Find RBH communities
#'
#' @description
#' This function will identify communities in the reciprocal best hit (RBH)
#' graph. It will iterate through resolutions and add the results to the
#' class. Additionally, a column will be added that signifies the resolution
#' with the best modularity.
#'
#' @param object The underlying class, see [bixverse::rbh_graph()].
#' @param resolution_params List. Parameters for the resolution search, see
#' [bixverse::params_graph_resolution()]. Contains:
#' \itemize{
#'  \item min_res - Float. Minimum resolution to test.
#'  \item max_res - Float. Maximum resolution to test.
#'  \item number_res - Integer. Number of resolutions to test between the
#'  `max_res` and `min_res.`
#' }
#' @param parallel Boolean. Shall the resolution search be in parallel.
#' @param max_workers Integer. Number of maximum cores to use. Defaults to half
#' of the identified cores (to a maximum of 8).
#' @param random_seed Integer. Random seed for reproducibility.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @return The class with added community detection results.
#'
#' @export
find_rbh_communities <- S7::new_generic(
  name = "find_rbh_communities",
  dispatch_args = "object",
  fun = function(
    object,
    resolution_params = params_graph_resolution(),
    max_workers = NULL,
    parallel = TRUE,
    random_seed = 42L,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)


#' @export
#'
#' @importFrom magrittr `%>%`
#' @import data.table
#'
#' @method find_rbh_communities rbh_graph
S7::method(find_rbh_communities, rbh_graph) <- function(
  object,
  resolution_params = params_graph_resolution(),
  max_workers = NULL,
  parallel = TRUE,
  random_seed = 42L,
  .verbose = TRUE
) {
  # Scope checks
  best_modularity <- modularity <- . <- NULL
  # Checks
  checkmate::assertClass(object, "bixverse::rbh_graph")
  assertGraphResParams(resolution_params)
  checkmate::qassert(parallel, "B1")
  checkmate::qassert(max_workers, c("I1", "0"))
  checkmate::qassert(.verbose, "B1")

  if (is.null(S7::prop(object, "rbh_graph"))) {
    warning("No RBH graph yet generated. Returning class as is.")
    return(object)
  }

  graph <- S7::prop(object, "rbh_graph")

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

  community_df_res <- data.table::rbindlist(community_df_res) %>%
    .[, best_modularity := modularity == max(modularity)]

  community_df_res[,
    c("origin_id", "module_id") := sub(
      "^([^_]+_[^_]+)_(.*)$",
      "\\1|\\2",
      community_df_res$node_name
    ) %>%
      strsplit(., "\\|") %>%
      transpose()
  ]

  community_df_res <- community_df_res[,
    c(
      "node_name",
      "origin_id",
      "module_id",
      "membership",
      "resolution",
      "modularity",
      "best_modularity"
    ),
    with = FALSE
  ]

  S7::prop(object, "final_results") <- community_df_res

  return(object)
}


# helpers ----------------------------------------------------------------------

## plots -----------------------------------------------------------------------

#' @export
#'
#' @import ggplot2
#'
#' @method plot_resolution_res rbh_graph
S7::method(plot_resolution_res, rbh_graph) <- function(
  object,
  print_head = TRUE,
  ...
) {
  checkmate::assertClass(object, "bixverse::rbh_graph")
  # Ignoring print_head for this class

  plot_df <- S7::prop(object, "final_results")
  if (is.null(plot_df)) {
    warning(
      paste(
        "No resolution results found. Did you run cor_module_check_res()?",
        "Returning NULL."
      )
    )
    return(NULL)
  }
  plot_df <- plot_df[, c("resolution", "modularity")] %>%
    unique()
  p <- ggplot2::ggplot(
    data = plot_df,
    mapping = aes(x = resolution, y = modularity)
  ) +
    ggplot2::geom_point(size = 3, shape = 21, alpha = .7) +
    ggplot2::xlab("Leiden cluster resolution") +
    ggplot2::ylab("Modularity") +
    ggplot2::theme_minimal() +
    ggplot2::ggtitle("Resolution vs. modularity")

  p
}
