# Methods ----

## Diffusion methods ----

#' Diffuse seed genes over network
#'
#' @description
#' This function takes a diffusion vector and leverages personalised page-rank diffusion
#' to identify influential nodes. These can be used subsequently for community detection
#' or check AUROC values given a set of genes.
#'
#' @usage diffuse_seed_genes(
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
#' @export
#'
#' @importFrom magrittr `%>%`
diffuse_seed_genes <- S7::new_generic("diffuse_seed_genes", "network_diffusions")

S7::method(diffuse_seed_genes, network_diffusions) <-
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
    seed_genes <- intersect(names(diffusion_vector), nodes_names)
    diff_vec <- rep(0, length(nodes_names)) %>% `names<-`(nodes_names)
    for (node in seed_genes) {
      diff_vec[node] <- diffusion_vector[node]
    }

    if (sum(diff_vec) == 0)
      error("No scores found to diffuse over the network. Please check the names and/or values.")

    ## Create the page-rank diffusion
    page_rank_score <- igraph::page_rank(S7::prop(network_diffusions, "graph"), personalized = diff_vec)

    ## Assign and return
    S7::prop(network_diffusions, "diffusion_res") <- page_rank_score$vector
    S7::prop(network_diffusions, "params")['diffusion_type'] <- 'single'
    S7::prop(network_diffusions, "params")[['seed_genes']] <- seed_genes
    network_diffusions
  }


#' Use tied diffusion for two sets of seed genes.
#'
#' @description
#' This function takes two sets of diffusion vector and leverages tied diffusion to identify an intersection
#' of influential nodes.
#'
#' @usage diffuse_seed_genes(
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
#' @export
#'
#' @importFrom magrittr `%>%`
tied_diffusion <- S7::new_generic("tied_diffusion", "network_diffusions")

S7::method(tied_diffusion, network_diffusions) <-
  function(network_diffusions,
           diffusion_vector.1,
           diffusion_vector.2,
           summarisation = c("max", "mean", "harmonic_sum"),
           score_aggregation = c("min", "max", "mean"),
           .verbose = TRUE) {
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
    seed_genes.1 <- intersect(names(diffusion_vector.1), nodes_names)
    seed_genes.2 <- intersect(names(diffusion_vector.2), nodes_names)
    diff_vec.1 <- diff_vec.2 <- rep(0, length(nodes_names)) %>% `names<-`(nodes_names)
    for (node in seed_genes.1) {
      diff_vec.1[node] <- diffusion_vector.1[node]
    }
    for (node in seed_genes.2) {
      diff_vec.2[node] <- diffusion_vector.2[node]
    }
    if ((sum(diff_vec.1) == 0) || (sum(diff_vec.1) == 0) )
      error("No scores found on first and/or second of the diffusion vectors. Please check the names and/or values.")

    ## First diffusion
    score_1 <- igraph::page_rank(S7::prop(network_diffusions, "graph"), personalized = diff_vec.1)$vector

    ## Second diffusion
    directed <- S7::prop(network_diffusions, "params")[['directed_graph']]
    score_2 <- if (directed) {
      if (.verbose)
        message("Directed graph found. Function will use transpose of adjacency for second diffusion.")
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
    S7::prop(network_diffusions, "params")[['seed_genes']] <- list('set_1' = seed_genes.1, 'set_2' = seed_genes.2)
    network_diffusions
  }

## Community detection ----

# Utils ----

## Others ----

calculate_diff_AUC <- S7::new_generic("calculate_diff_AUC", "network_diffusions")

S7::method(calculate_diff_AUC, network_diffusions) <-
  function(network_diffusions,
           hit_nodes,
           random_aucs = 1000L,
           auc_iters = 10000L,
           seed = 42L) {
    # Checks
    checkmate::assertClass(network_diffusions, "BIXverse::network_diffusions")
    checkmate::qassert(hit_nodes, "S+")
    checkmate::qassert(random_aucs, "I1")
    checkmate::qassert(auc_iters, "I1")
    checkmate::qassert(seed, "I1")
    # Body
    diffusion_score <- S7::prop(network_diffusions, "diffusion_res")
    if(length(diffusion_score) == 0) {
      warning("The diffusion score has length 0. Likely you did not run the diffusion methods. Returning NULL.")
      return(NULL)
    }
    pos.scores <- diffusion_score[hit_nodes]
    neg.scores <- diffusion_score[which(!names(scores) %in% hit_nodes)]
    auc <- rs_fast_auc(
      pos_scores = pos.scores,
      neg_scores = neg.scores,
      iters = auc_iters,
      seed = seed
    )
    random_aucs <- rs_create_random_aucs(
      score_vec = diffusion_score,
      size_pos = length(pos.scores),
      random_iters = random_aucs,
      auc_iters = auc_iters,
      seed = seed
    )

    z <- (auc - mean(random_aucs)) / sd(random_aucs)

    list(auc = auc, z = z)
  }


## Helpers ----

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
