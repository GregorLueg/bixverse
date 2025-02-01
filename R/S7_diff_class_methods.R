# Methods ----

## Diffusion methods ----

#' Diffuse seed genes over a network
#'
#' @description
#' This is the generic function for diffusing seed genes over a network.
#'
#' @export
diffuse_seed_nodes <- S7::new_generic("diffuse_seed_nodes", "network_diffusions")


#' @name diffuse_seed_nodes
#'
#' @description
#' This function takes a diffusion vector and leverages personalised page-rank diffusion
#' to identify influential nodes. These can be used subsequently for community detection
#' or check AUROC values given a set of genes.
#'
#' @usage diffuse_seed_nodes(
#'  network_diffusions,
#'  diffusion_vector,
#'  summarisation = c("max", "mean", "harmonic_sum")
#' )
#'
#' @param network_diffusions The underlying `network_diffusions` class.
#' @param diffusion_vector A named vector with values to use for the reset parameter in the
#' personalised page-rank diffusion. Names should represent node names of the graph.
#' @param summarisation If there are duplicated names in the `diffusion_vector` how to summarise
#' these.
#'
#' @return The class with added diffusion score based on a single set of seed genes. Additionally,
#' the seed genes are stored in the class.
#'
#' @export
#'
#' @importFrom magrittr `%>%`
#'
#' @method diffuse_seed_nodes network_diffusions
S7::method(diffuse_seed_nodes, network_diffusions) <-
  function(network_diffusions,
           diffusion_vector,
           summarisation = c("max", "mean", "harmonic_sum")) {
    # Checks
    checkmate::assertClass(network_diffusions, "BIXverse::network_diffusions")
    checkmate::assertNumeric(diffusion_vector)
    checkmate::assertNamed(diffusion_vector, .var.name = "diffusion_vector")
    checkmate::assertChoice(summarisation, c("max", "mean", "harmonic_sum"))

    # Body
    ## Create the diffusion vector
    diffusion_vector <- .summarise_scores(diffusion_vector, summarisation = summarisation)
    nodes_names <- igraph::V(S7::prop(network_diffusions, "graph"))$name
    seed_nodes <- intersect(names(diffusion_vector), nodes_names)
    diff_vec <- rep(0, length(nodes_names)) %>% `names<-`(nodes_names)
    for (node in seed_nodes) {
      diff_vec[node] <- diffusion_vector[node]
    }

    if (sum(diff_vec) == 0)
      error("No scores found to diffuse over the network. Please check the names and/or values.")

    ## Create the page-rank diffusion
    page_rank_score <- igraph::page_rank(S7::prop(network_diffusions, "graph"), personalized = diff_vec)

    ## Assign and return
    S7::prop(network_diffusions, "diffusion_res") <- page_rank_score$vector
    S7::prop(network_diffusions, "params")['diffusion_type'] <- 'single'
    S7::prop(network_diffusions, "params")[['seed_nodes']] <- seed_nodes
    network_diffusions
  }

#' Diffuse seed genes in a tied manner over a network
#'
#' @description
#' This is the generic function for diffusion two sets of seed genes over
#' a network.
#'
#' @export
tied_diffusion <- S7::new_generic("tied_diffusion", "network_diffusions")


#' @name tied_diffusion
#'
#' @description
#' This function takes two sets of diffusion vector and leverages tied diffusion to identify an intersection
#' of influential nodes.
#'
#' @usage tied_diffusion(
#'  network_diffusions,
#'  diffusion_vector.1,
#'  diffusion_vector.2,
#'  summarisation = c("max", "mean", "harmonic_sum"),
#'  score_aggregation = c("min", "max", "mean"),
#'  verbose = TRUE
#' )
#'
#' @param network_diffusions The underlying `network_diffusions` class.
#' @param diffusion_vector.1 The first named vector with values to use for the reset parameter in the
#' personalised page-rank diffusion. Names should represent node names of the graph.
#' @param diffusion_vector.2 The second named vector with values to use for the reset parameter in the
#' personalised page-rank diffusion. Names should represent node names of the graph.
#' @param summarisation If there are duplicated names in the `diffusion_vector` how to summarise
#' these.
#' @param score_aggregation How to summarise the tied scores.
#' @param verbose Controls verbosity of the function.
#'
#' @return The class with added diffusion score based on a two sets of seed genes. Additionally,
#' the seed genes are stored in the class.
#'
#' @export
#'
#' @importFrom magrittr `%>%`
#'
#' @method tied_diffusion network_diffusions
S7::method(tied_diffusion, network_diffusions) <-
  function(network_diffusions,
           diffusion_vector.1,
           diffusion_vector.2,
           summarisation = c("max", "mean", "harmonic_sum"),
           score_aggregation = c("min", "max", "mean"),
           .verbose = FALSE) {
    # Checks
    checkmate::assertClass(network_diffusions, "BIXverse::network_diffusions")
    checkmate::assertNumeric(diffusion_vector.1)
    checkmate::assertNamed(diffusion_vector.1, .var.name = "diffusion_vector.1")
    checkmate::assertNumeric(diffusion_vector.2)
    checkmate::assertNamed(diffusion_vector.2, .var.name = "diffusion_vector.2")
    checkmate::assertChoice(summarisation, c("max", "mean", "harmonic_sum"))
    checkmate::assertChoice(score_aggregation, c("min", "max", "mean"))
    checkmate::qassert(.verbose, "B1")

    # Body
    ## Create the diffusion vectors
    diffusion_vector.1 <- .summarise_scores(diffusion_vector.1, summarisation = summarisation)
    diffusion_vector.2 <- .summarise_scores(diffusion_vector.2, summarisation = summarisation)
    nodes_names <- igraph::V(S7::prop(network_diffusions, "graph"))$name
    seed_nodes.1 <- intersect(names(diffusion_vector.1), nodes_names)
    seed_nodes.2 <- intersect(names(diffusion_vector.2), nodes_names)
    diff_vec.1 <- diff_vec.2 <- rep(0, length(nodes_names)) %>% `names<-`(nodes_names)
    for (node in seed_nodes.1) {
      diff_vec.1[node] <- diffusion_vector.1[node]
    }
    for (node in seed_nodes.2) {
      diff_vec.2[node] <- diffusion_vector.2[node]
    }
    if ((sum(diff_vec.1) == 0) || (sum(diff_vec.1) == 0))
      error(
        "No scores found on first and/or second of the diffusion vectors. Please check the names and/or values."
      )

    ## First diffusion
    score_1 <- igraph::page_rank(S7::prop(network_diffusions, "graph"), personalized = diff_vec.1)$vector

    ## Second diffusion
    directed <- S7::prop(network_diffusions, "params")[['directed_graph']]
    score_2 <- if (directed) {
      if (.verbose)
        message(
          "Directed graph found. Function will use transpose of adjacency for second diffusion."
        )
      adj <- igraph::as_adjacency_matrix(S7::prop(network_diffusions, "graph"))
      adj_T <- Matrix::t(adj)
      igraph_obj_t <- igraph::graph_from_adjacency_matrix(adj_T)
      igraph::page_rank(igraph_obj_t, personalized = diff_vec.2)$vector
    } else {
      if (.verbose)
        message("Undirected graph found. Using graph as is for second diffusion.")
      igraph::page_rank(S7::prop(network_diffusions, "graph"), personalized = diff_vec.2)$vector
    }

    ## Summarise the scores
    final_tiedie_diffusion <- switch(
      score_aggregation,
      "min" = pmin(score_1, score_2),
      "max" = pmax(score_1, score_2),
      rowMeans(cbind(score_1, score_2))
    )

    ## Assign and return
    S7::prop(network_diffusions, "diffusion_res") <- final_tiedie_diffusion
    S7::prop(network_diffusions, "params")['diffusion_type'] <- 'tied'
    S7::prop(network_diffusions, "params")[['seed_nodes']] <- list('set_1' = seed_nodes.1, 'set_2' = seed_nodes.2)

    network_diffusions
  }

## Community detection ----

#' Identify privileged communities based on a given diffusion vector
#'
#' @description
#' This is the generic function for diffusing seed genes over a network.
#'
#' @export
community_detection <- S7::new_generic("community_detection", "network_diffusions")

#' @name community_detection
#'
#' @description Detects privileged communities after a diffusion based on seed nodes.
#'
#' @usage community_detection(
#'  network_diffusions,
#'  diffusion_threshold,
#'  max_nodes = 300L,
#'  min_nodes = 10L,
#'  min_seed_nodes = 2L,
#'  intial_res = 0.5,
#'  seed = 42L,
#'  .verbose = F,
#'  .max_iters = 100L
#' )
#'
#' @param network_diffusions The underlying `network_diffusions` class.
#' @param diffusion_threshold How much of the network to keep based on the diffusion
#' values. 0.25 for example would keep the 25% nodes with the highest scores.
#' @param max_nodes Number of max genes per community. Communities that are larger
#' than that will be further subclustered.
#' @param min_nodes Number of minimum nodes per community. Smaller communities will be removed.
#' @param min_seed_nodes Number of minimum seed nodes per community to be kept. In the case of
#' a tied diffusion, this number must be reached by both initial seed gene vectors.
#' @param intial_res Initial resolution parameter for the Leiden community detection.
#' @param seed Random seed.
#' @param .verbose Controls the verbosity of the function.
#' @param .max_iters Controls how many iterations shall be tried for the subclustering. To note,
#' in each iteration of the subclustering, the resolution parameter is increased by 0.05, to identify
#' more granular communities.
#'
#' @return The class with added diffusion community detection results (if any could be identified
#' with the provided parameters).
#'
#' @export
#'
#' @import data.table
#' @importFrom magrittr `%>%`
#'
#' @method community_detection network_diffusions
S7::method(community_detection, network_diffusions) <- function(network_diffusions,
                                                                diffusion_threshold,
                                                                max_nodes = 300L,
                                                                min_nodes = 10L,
                                                                min_seed_nodes = 2L,
                                                                intial_res = 0.5,
                                                                seed = 42L,
                                                                .verbose = F,
                                                                .max_iters = 100L) {
  # Checks
  checkmate::assertClass(network_diffusions, "BIXverse::network_diffusions")
  checkmate::qassert(diffusion_threshold, "R1[0,1]")
  checkmate::qassert(max_nodes, "I1")
  checkmate::qassert(min_nodes, "I1")
  checkmate::qassert(min_seed_nodes, "I1")
  checkmate::qassert(intial_res, "R1")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, "B1")
  checkmate::qassert(.max_iters, "I1")

  # Body
  ## Reduce the graph
  diffusion_score <- S7::prop(network_diffusions, "diffusion_res")
  # Early return
  if (length(diffusion_score) == 0) {
    warning(
      "The diffusion score has length 0. Likely you did not run the diffusion methods. Returning class as is."
    )
    return(network_diffusions)
  }
  nodes_to_include <- diffusion_score %>%
    sort(decreasing = T) %>%
    .[1:ceiling(diffusion_threshold * length(diffusion_score))]

  red_graph <- igraph::subgraph(S7::prop(network_diffusions, "graph"),
                                names(nodes_to_include))

  ## First clustering
  set.seed(seed)

  if (S7::prop(network_diffusions, "params")$directed_graph) {
    red_graph = igraph::as_undirected(red_graph, mode = 'each')
  }

  first_clusters <- igraph::cluster_leiden(red_graph, resolution = intial_res)

  clusters_df <- data.table::data.table(nodeID = first_clusters$names, clusterID = first_clusters$membership)

  node_frequency <- clusters_df[, .N, .(clusterID)]

  ## Subclustering
  clusters_with_too_many_nodes <- node_frequency[N > max_nodes, clusterID]
  final_clusters <- clusters_df[!clusterID %in% clusters_with_too_many_nodes]

  for (i in seq_along(clusters_with_too_many_nodes)) {
    cluster_i <- clusters_with_too_many_nodes[i]
    nodes_in_cluster <- clusters_df[clusterID == cluster_i, nodeID]
    finalised_clusters <- data.table()
    # Loop through, until all clusters are below the minimum genes or max iterations is hit
    l <- 1
    while (length(nodes_in_cluster) != 0) {
      set.seed(seed + l)
      if (.verbose)
        message("Cluster ", i, " gets subclustered. Iter: ", l)
      red_graph_l <- igraph::subgraph(red_graph,
                                      data.table::chmatch(nodes_in_cluster, igraph::V(red_graph)$name))

      clusters_red <- igraph::cluster_leiden(red_graph_l, resolution = intial_res + l * 0.05)

      subclusters <- data.table(nodeID = clusters_red$names,
                                clusterID = clusters_red$membership)
      subclusters_frequency <- subclusters[, .N, .(clusterID)]
      clusters_small_enough <- subclusters_frequency[N <= max_nodes, clusterID]

      good_clusters <- subclusters[clusterID %in% clusters_small_enough] %>%
        dplyr::mutate(clusterID = paste0(i, paste(rep("sub", l), collapse = ""), clusterID))

      finalised_clusters <- rbind(finalised_clusters, good_clusters)

      l <- l + 1
      if (l == .max_iters)
        break

      nodes_in_cluster <- setdiff(nodes_in_cluster, good_clusters$nodeID)
    }

    final_clusters <- rbind(final_clusters, finalised_clusters)
  }

  ## Add the seed node information based on diffusion type
  diffusion_type <- S7::prop(network_diffusions, "params")$diffusion_type

  final_node_frequency <- if (diffusion_type == 'single') {
    seed_nodes <- S7::prop(network_diffusions, "params")$seed_nodes

    final_clusters[, .(
      cluster_size = length(nodeID),
      seed_nodes_no = sum(nodeID %in% seed_nodes)
    ), .(clusterID)]
  } else {
    seed_nodes_set_1 <- S7::prop(network_diffusions, "params")$seed_nodes$set_1
    seed_nodes_set_2 <- S7::prop(network_diffusions, "params")$seed_nodes$set_2

    final_clusters[, .(
      cluster_size = length(nodeID),
      seed_nodes_1 = sum(nodeID %in% seed_nodes_set_1),
      seed_nodes_2 = sum(nodeID %in% seed_nodes_set_2)
    ), .(clusterID)]
  }

  ## Finalise the clusters
  clusters_to_take <- if (diffusion_type == 'single') {
    final_node_frequency[cluster_size >= min_nodes &
                           seed_nodes_no >= min_seed_nodes, clusterID]
  } else {
    final_node_frequency[cluster_size >= min_nodes &
                           seed_nodes_1 >= min_seed_nodes &
                           seed_nodes_2 >= min_seed_nodes, clusterID]
  }

  # Early return
  if (length(clusters_to_take) == 0) {
    warning("No communities found with the given parameters. Returning class as is.")
    return(network_diffusions)
  }

  finalised_clusters.clean <- final_clusters[clusterID %in% clusters_to_take]

  ks_vals = vector(mode = "numeric", length = length(clusters_to_take))

  for (i in seq_along(clusters_to_take)) {
    cluster <- clusters_to_take[i]
    cluster_nodes <- finalised_clusters.clean[clusterID == cluster, nodeID]
    ks <- suppressWarnings(ks.test(diffusion_score[cluster_nodes], diffusion_score[which(!names(diffusion_score) %in% cluster_nodes)], alternative = "less"))
    ks_vals[i] <- ks$p.value
  }

  ks_val_df <- data.table(clusterID = clusters_to_take, ks_pval = ks_vals)

  final_result <- purrr::reduce(list(finalised_clusters.clean, ks_val_df, final_node_frequency),
                                merge,
                                by = 'clusterID') %>%
    .[, `:=`(diffusion_score = diffusion_score[nodeID],
             clusterID = as.character(clusterID))]


  cluster_name_prettifier <- setNames(paste("cluster", 1:length(unique(
    final_result$clusterID
  )), sep = "_"),
  unique(final_result$clusterID))

  final_result[, clusterID := cluster_name_prettifier[clusterID]]

  ## Assign and return
  S7::prop(network_diffusions, "community_res") <- final_result
  S7::prop(network_diffusions, "params")[['community_params']] <- list(
    diffusion_threshold = diffusion_threshold,
    max_nodes = max_nodes,
    min_nodes = min_seed_nodes,
    min_seed_nodes = min_seed_nodes
  )

  return(network_diffusions)
}

# Utils ----

## Statistics ----

#' Calculate the AUROC for a diffusion score
#'
#' @description
#' This is the generic function for calculating AUC/AUROCs for a set of genes based
#' on the diffusion vector of either single or tied diffusion.
#'
#' @export
calculate_diffusion_AUC <- S7::new_generic("calculate_diffusion_AUC", "network_diffusions")


#' @name calculate_diffusion_AUC
#'
#' @description
#' This functions can take a given `network_diffusions` class and calculates an AUC and generates
#' a Z-score based on random permutation of `random_aucs` for test for statistical significance if
#' desired.
#'
#' @usage calculate_diffusion_AUC(
#'  network_diffusions,
#'  hit_nodes,
#'  auc_iters = 10000L,
#'  random_aucs = 1000L,
#'  permutation_test = FALSE,
#'  seed = 42L
#' )
#'
#' @param network_diffusions The underlying `network_diffusions` class.
#' @param hit_nodes Which nodes would consist a hit to calculate the AUROC against.
#' @param auc_iters How many iterations to run to approximate the AUROC.
#' @param random_aucs How many random AUROCs to calculate to estimate the Z-score
#' @param seed Seed for reproducibility purposes.
#'
#' @return List with AUC and Z-score as the two named elements if permutations_test set
#' to TRUE; otherwise just the AUC.
#'
#' @export
#'
#' @importFrom magrittr `%>%`
#'
#' @method tied_diffusion network_diffusions
S7::method(calculate_diffusion_AUC, network_diffusions) <-
  function(network_diffusions,
           hit_nodes,
           auc_iters = 10000L,
           random_aucs = 1000L,
           permutation_test = FALSE,
           seed = 42L) {
    # Checks
    checkmate::assertClass(network_diffusions, "BIXverse::network_diffusions")
    checkmate::qassert(hit_nodes, "S+")
    checkmate::qassert(auc_iters, "I1")
    checkmate::qassert(seed, "I1")
    checkmate::qassert(permutation_test, "B1")
    if (permutation_test)
      checkmate::qassert(random_aucs, "I1")
    # Body
    diffusion_score <- S7::prop(network_diffusions, "diffusion_res")
    if (length(diffusion_score) == 0) {
      warning(
        "The diffusion score has length 0. Likely you did not run the diffusion methods. Returning NULL."
      )
      return(NULL)
    }
    pos.scores <- diffusion_score[hit_nodes]
    neg.scores <- diffusion_score[which(!names(diffusion_score) %in% hit_nodes)]
    auc <- rs_fast_auc(
      pos_scores = pos.scores,
      neg_scores = neg.scores,
      iters = auc_iters,
      seed = seed
    )
    if (permutation_test) {
      random_aucs <- rs_create_random_aucs(
        score_vec = diffusion_score,
        size_pos = length(pos.scores),
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

# Other Helpers ----

#' Summarise gene scores if they are duplicates.
#'
#' @param x Named numeric.
#' @param summarisation String. Which summary function to use.
#'
#' @return Named numeric.
#'
#' @export
#'
#' @importFrom magrittr `%$%`
.summarise_scores <- function(x,
                              summarisation = c("max", "mean", "harmonic_sum")) {
  # Checks
  checkmate::assertNumeric(x)
  checkmate::assertNamed(x, .var.name = "x")
  checkmate::assertChoice(summarisation, c("max", "mean", "harmonic_sum"))
  # Body
  dt <- data.table::data.table(node_name = names(x), value = x)
  summary_fun <- if (summarisation == "mean") {
    rlang::expr(mean(value))
  } else if (summarisation == "max") {
    rlang::expr(max(value))
  } else {
    rlang::expr(BIXverse::OT_harmonic_score(value))
  }
  res <-
    rlang::eval_tidy(rlang::quo(dt[, .(value = !!summary_fun), .(node_name)])) %$%
    setNames(value, node_name)
  res
}
