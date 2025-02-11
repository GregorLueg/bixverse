# S7 objects ----

## network_diffusions ----

### class ----

network_diffusions <- S7::new_class(
  # Name
  name = "network_diffusions",
  parent = bixverse_generic_class,
  # Properties, i.e., slots
  properties = list(
    graph = S7::class_any,
    diffusion_res = S7::class_numeric,
    final_results = S7::class_data.frame,
    params = S7::class_list
  ),

  #' Network diffusion class
  #'
  #' @description
  #' This class helps to do diffusion of seed nodes in a single or tied version
  #' over a network, measure the ability of these diffusion vectors to recall
  #' against a gold standard set of nodes and do community detection within the
  #' subset of the network that received the the most heat from the initial seed
  #' genes.
  #'
  #' @param edge_data_frame data.table that contains the edge information. It is
  #' expected to have the columns 'from' and 'to'.
  #' @param weighted Boolean. Is the graph weighted. If set to TRUE, the
  #' `edge_data_frame` needs to have a weight column.
  #' @param directed Boolean. Shall the graph be stored as directed.
  #'
  #' @return Returns the `network_diffusions` class for further operations.
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
    params <- list(
      "directed_graph" = igraph::is_directed(graph),
      "weighted_graph" = igraph::is_weighted(graph)
    )

    # Finalise object
    S7::new_object(
      S7::S7_object(),
      graph = graph,
      diffusion_res = vector(mode = "numeric"),
      final_results = data.table(),
      params = params
    )
  }
)

### methods ----



### print/show ----

S7::method(print, network_diffusions) <- function(x, ...) {
  x <- length(S7::prop(x, "diffusion_res"))

  cat(sprintf("Diffusion vector has length %i", x))
}


## rbh_graphs ----

rbh_graph <- S7::new_class(
  # Name
  name = "rbh_graph",
  parent = bixverse_generic_class,
  # Properties, i.e., slots
  properties = list(
    module_data = S7::class_list,
    rbh_graph = S7::class_any,
    rbh_edge_df = S7::class_data.frame,
    final_results = S7::class_data.frame,
    params = S7::class_list
  ),

  #' Reciprocal best hit graph
  #'
  #' @description
  #' This class can be used to generate reciprocal best hit graphs between
  #' gene modules from different origins.
  #'
  #' @param module_results data.table with the all of the gene modules for which you
  #' wish to generate the RBH graph.
  #' @param dataset_col The column (name) which indicates from which data set/method
  #' the gene module was derived.
  #' @param module_col The column (name) which stores the names of the modules.
  #' @param value_col The column (name) which stores the genes that are part of the modules.
  #'
  #' @return Returns the `rbh_graph` class for further operations.
  #'
  #' @export
  constructor = function(module_results,
                         dataset_col,
                         module_col,
                         value_col) {
    checkmate::assertDataTable(module_results)
    checkmate::qassert(dataset_col, "S1")
    checkmate::qassert(module_col, "S1")
    checkmate::qassert(value_col, "S1")
    checkmate::assertNames(names(module_results),
      must.include = c(dataset_col, module_col, value_col)
    )
    # Function body
    list_of_list <- split(
      module_results %>% dplyr::select(!!module_col, !!value_col),
      module_results[, ..dataset_col]
    ) %>%
      purrr::map(., ~ {
        df <- .
        split(unlist(df[, ..value_col]), unlist(df[, ..module_col]))
      })

    # Finalise object
    S7::new_object(
      S7::S7_object(),
      module_data = list_of_list,
      rbh_edge_df = data.table(),
      rbh_graph = NULL,
      final_results = data.table(),
      params = list(
        no_compared_modules = nrow(module_results)
      )
    )
  }
)

## methods ----

### print/show ----



### getters ----
