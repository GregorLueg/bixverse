# S7 ----

network_diffusions <- S7::new_class(
  # Name
  name = "network_diffusions",
  # Properties, i.e., slots
  properties = list(
    graph = S7::class_any,
    diffusion_res = S7::class_numeric,
    community_res = S7::class_data.frame,
    params = S7::class_list
  ),

  #' Network diffusion class
  #'
  #' @description
  #' ...
  #'
  #' @param edge_data_frame data.table with the columns
  #' @param weighted Boolean. Is the graph weighted.
  #' @param directed Boolean. Shall the graph be stored as directed.
  #'
  #' @return Returns the S7 object for further operations.
  #'
  #' @export
  constructor = function(edge_data_frame, weighted, directed) {
    # Checks
    needed_cols <- if (weighted) {
      c("from", "to", "weight")
    } else {
      c("from", "to")
    }
    checkmate::assertDataTable(edge_data_frame)
    checkmate::assertNames(names(edge_data_frame), must.include = needed_cols)
    checkmate::qassert(weighted, "B1")
    checkmate::qassert(directed, "B1")
    # Function body
    graph <- igraph::graph_from_data_frame(edge_data_frame, directed = directed)
    params <-list(
      "directed_graph" = igraph::is_directed(graph),
      "weighted_graph" = igraph::is_weighted(graph)
    )

    # Finalise object
    S7::new_object(
      S7::S7_object(),
      graph = graph,
      diffusion_res = vector(mode = 'numeric'),
      community_res = data.table(),
      params = params
    )
  }
)
