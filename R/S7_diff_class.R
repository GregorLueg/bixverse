# S7 object ----

## Class ----

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
  #' This class helps to do diffusion of seed nodes in a single or tied version over a
  #' network, measure the ability of these diffusion vectors to recall against a gold standard
  #' set of nodes and do community detection within the subset of the network that received the
  #' the most heat from the initial seed genes.
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

## Property access methods ----

#' Get the parameters that were used.
#'
#' @description
#' This method accesses the params slot and can return R lists or JSON strings.
#'
#' @export
get_params <- S7::new_generic("get_params", "network_diffusions")

#' @name get_params
#'
#' @description Extracts params from the `network_diffusion` class and has options
#' to return (pretty) JSONs
#'
#' @usage get_params(
#'  network_diffusions,
#'  to_json = FALSE,
#'  pretty_json = FALSE
#' )
#'
#' @param network_diffusions The underlying `network_diffusions` class.
#' @param to_json Shall the params be returned as a JSON string.
#' @param pretty_json Shall the params be returned as a pretty JSON string.
#'
#' @return Depending on parameters either the R list or a (pretty) JSON string.
#'
#' @method get_params network_diffusions
S7::method(get_params, network_diffusions) <-
  function(network_diffusions,
           to_json = FALSE,
           pretty_json = FALSE) {
    # Checks
    checkmate::assertClass(network_diffusions, "BIXverse::network_diffusions")
    checkmate::qassert(to_json, "B1")
    checkmate::qassert(pretty_json, "B1")

    # Body
    to_ret <- S7::prop(network_diffusions, "params")
    if (to_json) {
      to_ret <- jsonlite::toJSON(to_ret)
    }
    if (to_json &&
        pretty_json) {
      to_ret <- jsonlite::prettify(to_ret)
    }

    return(to_ret)
  }


#' Get the diffusion results
#'
#' @description
#' This method returns the community detection results from the class
#'
#' @export
get_results <- S7::new_generic("get_results", "network_diffusions")

#' @name get_results
#'
#' @description Get the community detection results from the class
#'
#' @usage get_results(network_diffusions)
#'
#' @param network_diffusions The underlying `network_diffusions` class.
#'
#' @return Returns the community detection results if any can be found.
#'
#' @method get_results network_diffusions
S7::method(get_results, network_diffusions) <-
  function(network_diffusions) {
    # Checks
    checkmate::assertClass(network_diffusions, "BIXverse::network_diffusions")

    # Return
    return(S7::prop(network_diffusions, "community_res"))
  }
