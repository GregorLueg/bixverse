# constrained page rank --------------------------------------------------------

#' Constrained personalised page rank
#'
#' @description
#' This function implements a constrained personalised page-rank. You can define
#' `sink_nodes` or `sink_edges`. In the case of the former, the vist of a
#' sink_node (via a type attribute in the igraph) will automatically cause the
#' surfer to reset. In the case of the latter, the traversal of a sink_edge is
#' allowed, however, subsequently, the surfer will be reseted. This function can
#' be useful for deriving PPR profiles under constraints in heterogenous graphs
#' and has been inspired by Ruiz, et al.
#'
#' @param graph igraph. This one needs to be directed and weighted and have
#' the node attribute `type` defining the node type and the edge attribute
#' `type` defining the edge type.
#' @param personalisation_vector Numerical vector. The personalisation vector
#' for the reset of the personalised page rank. The function will automatically
#' normalise it.
#' @param sink_nodes Optional String vector. The node types that should force
#' a reset.
#' @param sink_edges Optional String vector. The edge types after which there
#' should be a forced reset.
#'
#' @returns The constrained personalised page rank
#'
#' @export
#'
#' @references Ruiz, et al., Nat Commun, 2021
constrained_page_rank <- function(
  graph,
  personalisation_vector,
  sink_nodes = NULL,
  sink_edges = NULL
) {
  # checks
  checkmate::assertClass(graph, "igraph")
  checkmate::assertTRUE(igraph::is_directed(graph))
  checkmate::assertTRUE(igraph::is_weighted(graph))
  checkmate::assertTRUE(!is.null(igraph::V(graph)$type))
  checkmate::assertTRUE(!is.null(igraph::E(graph)$type))
  checkmate::assertNumeric(
    personalisation_vector,
    len = length(igraph::V(graph))
  )
  checkmate::qassert(sink_nodes, c("S+", "0"))
  checkmate::qassert(sink_edges, c("S+", "0"))

  # function body
  edge_list <- igraph::as_edgelist(graph, names = TRUE)
  personalisation_vector <- personalisation_vector / sum(personalisation_vector)

  res <- rs_constrained_page_rank(
    node_names = igraph::vertex.attributes(graph)$name,
    node_types = igraph::vertex.attributes(graph)$type,
    from = edge_list[, 1],
    to = edge_list[, 2],
    weights = igraph::edge.attributes(graph)$weight,
    edge_type = igraph::edge.attributes(graph)$type,
    personalised = personalisation_vector,
    sink_nodes = sink_nodes,
    sink_edges = sink_edges
  )
  names(res) <- igraph::vertex.attributes(graph)$name

  return(res)
}

#' Constrained personalised page rank over a list
#'
#' @description
#' This function implements a personalised constrained page-rank over a list of
#' personalisation vectors for the same network. You can define `sink_nodes`
#' or `sink_edges`. In the case of the former, the vist of a sink_node (via a
#' type attribute in the igraph) will automatically cause the surfer to reset.
#' In the case of the latter, the traversal of a sink_edge is allowed, however,
#' subsequently, the surfer will be reseted. This function can be useful for
#' deriving PPR profiles under constraints in heterogenous graphs and has been
#' inspired by Ruiz, et al.
#'
#' @param graph igraph. This one needs to be directed and weighted and have
#' the node attribute `type` defining the node type and the edge attribute
#' `type` defining the edge type.
#' @param personalisation_list A list of numerical vectors to use for the
#' personalisation.
#' @param sink_nodes Optional String vector. The node types that should force
#' a reset.
#' @param sink_edges Optional String vector. The edge types after which there
#' should be a forced reset.
#'
#' @returns A list with the constrained personalised page rank values.
#'
#' @export
#'
#' @references Ruiz, et al., Nat Commun, 2021
constrained_page_rank_ls <- function(
  graph,
  personalisation_list,
  sink_nodes = NULL,
  sink_edges = NULL
) {
  # checks
  checkmate::assertClass(graph, "igraph")
  checkmate::assertTRUE(igraph::is_directed(graph))
  checkmate::assertTRUE(igraph::is_weighted(graph))
  checkmate::assertTRUE(!is.null(igraph::V(graph)$type))
  checkmate::assertTRUE(!is.null(igraph::E(graph)$type))
  checkmate::assertList(personalisation_list, types = "numeric")
  checkmate::assertTRUE(
    all(purrr::map_dbl(personalisation_list, length) == length(graph))
  )
  checkmate::qassert(sink_nodes, c("S+", "0"))
  checkmate::qassert(sink_edges, c("S+", "0"))

  edge_list <- igraph::as_edgelist(graph, names = TRUE)
  # Ensure everything sums to 1
  personalisation_list <- purrr::map(personalisation_list, \(pers_vec) {
    pers_vec / sum(pers_vec)
  })

  results <- rs_constrained_page_rank_list(
    personalisation_list = personalisation_list,
    node_names = igraph::vertex.attributes(graph)$name,
    node_types = igraph::vertex.attributes(graph)$type,
    from = edge_list[, 1],
    to = edge_list[, 2],
    weights = igraph::edge.attributes(graph)$weight,
    edge_type = igraph::edge.attributes(graph)$type,
    sink_nodes = sink_nodes,
    sink_edges = sink_edges
  )
  # add all the names
  names(results) <- names(personalisation_list)
  results <- purrr::map(
    results,
    ~ {
      names(.x) <- igraph::vertex.attributes(graph)$name
      .x
    }
  )

  return(results)
}

## helpers ---------------------------------------------------------------------

#' Helper function to create personalisation vectors
#'
#' @param graph igraph. The graph for which to produce the personalisation
#' vector.
#' @param node_weights Named numeric. The names represent the nodes and the
#' values the strength of the reset.
#'
#' @returns The personalisation vector for subsequent usage in page-rank
#'
#' @export
#'
#' @importFrom magrittr `%>%`
generate_personalisation_vec <- function(graph, node_weights) {
  # checks
  checkmate::assertClass(graph, "igraph")
  checkmate::qassert(node_weights, "N1")
  checkmate::assertNamed(node_weights)

  # function body
  intersecting_nodes <- intersect(names(node_weights), igraph::V(graph)$name)
  diffusion_vec <- rep(0, length(graph)) %>% `names<-`(igraph::V(graph)$name)
  diffusion_vec[intersecting_nodes] <- node_weights[intersecting_nodes]
  if (sum(diffusion_vec) == 0) {
    stop(paste(
      "The function could not identify any intersecting nodes.",
      "Check the names please."
    ))
  }
  diffusion_vec < diffusion_vec / sum(diffusion_vec)

  return(diffusion_vec)
}
