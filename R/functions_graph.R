# constrained page rank --------------------------------------------------------

#' Constrained personalised page rank
#'
#' @description
#' This function implements a constrained page-rank. You can define `sink_nodes`
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
    node_names = igraph::V(graph)$name,
    node_types = igraph::V(graph)$type,
    from = edge_list[, 1],
    to = edge_list[, 2],
    weights = igraph::E(graph)$weight,
    edge_type = igraph::E(graph)$type,
    personalised = personalisation_vector,
    sink_nodes = sink_nodes,
    sink_edges = sink_edges
  )

  res
}
