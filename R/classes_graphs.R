# network_diffusions -----------------------------------------------------------

## class -----------------------------------------------------------------------

#' Network diffusion class
#'
#' @description
#' This class helps to do diffusion of seed nodes in a single or tied version
#' over a network, measure the ability of these diffusion vectors to recall
#' against a gold standard set of nodes and do community detection within the
#' subset of the network that received the the most heat from the initial seed
#' genes.
#'
#' @section Properties:
#' \describe{
#'   \item{graph}{igraph. The underlying graph.}
#'   \item{diffusion_res}{Numeric vector. Contains contains the single or tied
#'   diffusion results.}
#'   \item{params}{A (nested) list that will store all the parameters of the
#'   applied function.}
#'   \item{final_results}{data.table. Contains final results.}
#' }
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
network_diffusions <- S7::new_class(
  # Names, parents
  name = "network_diffusions",
  parent = bixverse_base_class,

  # Properties, i.e., slots
  properties = list(
    graph = S7::class_any,
    diffusion_res = S7::class_numeric,
    diffusion_perm = S7::class_numeric,
    final_results = S7::class_data.frame,
    params = S7::class_list
  ),
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
    graph <- if (weighted) {
      igraph::graph_from_data_frame(edge_data_frame, directed = directed)
    } else {
      igraph::graph_from_data_frame(
        edge_data_frame[, c("from", "to")],
        directed = directed
      )
    }

    if (!weighted) {
      params <- list(
        "directed_graph" = igraph::is_directed(graph),
        "weighted_graph" = igraph::is_weighted(graph)
      )
    }

    # Finalise object
    S7::new_object(
      S7::S7_object(),
      graph = graph,
      diffusion_res = vector(mode = "numeric"),
      diffusion_perm = vector(mode = "numeric"),
      final_results = data.table(),
      params = params
    )
  }
)

## methods ---------------------------------------------------------------------

### getters --------------------------------------------------------------------

#' Get the diffusion vector
#'
#' @description Returns the diffusion vector if you ran [bixverse::tied_diffusion()]
#' or [bixverse::diffuse_seed_nodes()].
#'
#' @param object The underlying class [bixverse::network_diffusions()].
#'
#' @return The diffusion vector if found. If you did not run either diffusion
#' functions, it will return `NULL` and a warning.
#'
#' @export
get_diffusion_vector <- S7::new_generic(
  name = "get_diffusion_vector",
  dispatch_args = "object",
  fun = function(object) {
    S7::S7_dispatch()
  }
)

#' @export
#'
#' @importFrom magrittr `%>%`
#'
#' @method get_diffusion_vector network_diffusions
S7::method(get_diffusion_vector, network_diffusions) <- function(object) {
  # Checks
  checkmate::assertClass(object, "bixverse::network_diffusions")
  # Get the data
  diffusion_vec <- S7::prop(object, "diffusion_res")
  if (is.null(diffusion_vec)) {
    warning("No diffusion results found. Returning NULL.")
  }
  diffusion_vec
}

#' Get the diffusion permutations
#'
#' @description Returns the diffusion Z-scores if you ran
#' [bixverse::permute_seed_nodes()].
#'
#' @param object The underlying class [bixverse::network_diffusions()].
#'
#' @return The diffusion Z scores if found. Otherwise `NULL`.
#'
#' @export
get_diffusion_perms <- S7::new_generic(
  name = "get_diffusion_perms",
  dispatch_args = "object",
  fun = function(object) {
    S7::S7_dispatch()
  }
)

#' @export
#'
#' @importFrom magrittr `%>%`
#'
#' @method get_diffusion_perms network_diffusions
S7::method(get_diffusion_perms, network_diffusions) <- function(object) {
  # Checks
  checkmate::assertClass(object, "bixverse::network_diffusions")
  # Get the data
  diffusion_perm <- S7::prop(object, "diffusion_perm")
  if (is.null(diffusion_perm)) {
    warning("No diffusion permutations found. Returning NULL.")
  }
  diffusion_perm
}

# rbh_graphs -------------------------------------------------------------------

## class -----------------------------------------------------------------------

#' Reciprocal best hit graph
#'
#' @description
#' This class can be used to generate reciprocal best hit graphs between
#' gene modules from different origins.
#'
#' @section Properties:
#' \describe{
#'   \item{module_data}{Nested list of the modules generated by different
#'   methods.}
#'   \item{rbh_graph}{igraph. The reciprocal best hit graph.}
#'   \item{rbh_edge_df}{data.table. Contains the RBH graph data as an edge df.}
#'   \item{params}{A (nested) list that will store all the parameters of the
#'   applied function.}
#'   \item{final_results}{data.table. Contains final results.}
#' }
#'
#' @param module_results data.table or list. If data.table (you set the class
#' to `"set"`, i.e., set similarities will be used) it cointains all of
#' the gene modules for which you wish to generate the RBH graph and you need
#' to have the three column names. If list (you set the class to `"cor"``,
#' i.e., correlations between for example gene loadings will be used), the
#' list contains all of the matrices (must have col and row names).
#' @param rbh_type String. One of `c("cor", "set")`.
#' @param dataset_col The column (name) which indicates from which data set/
#' method the gene module was derived. Only needed if you want to use set
#' similarities.
#' @param module_col The column (name) which stores the names of the modules.
#' Only needed if you want to use set similarities.
#' @param value_col The column (name) which stores the genes that are part of
#' the modules. Only needed if you want to use set similarities.
#'
#' @return Returns the `rbh_graph` class for further operations.
#'
#' @export
rbh_graph <- S7::new_class(
  # Names, parents
  name = "rbh_graph",
  parent = bixverse_base_class,

  # Properties, i.e., slots
  properties = list(
    module_data = S7::class_list,
    rbh_graph = S7::class_any,
    rbh_edge_df = S7::class_data.frame,
    final_results = S7::class_data.frame,
    params = S7::class_list
  ),

  constructor = function(
    module_results,
    rbh_type = c("cor", "set"),
    dataset_col = NULL,
    module_col = NULL,
    value_col = NULL
  ) {
    # checks
    checkmate::assertChoice(rbh_type, c("cor", "set"))
    if (rbh_type == "set") {
      checkmate::assertDataTable(module_results)
      checkmate::qassert(dataset_col, "S1")
      checkmate::qassert(module_col, "S1")
      checkmate::qassert(value_col, "S1")
      checkmate::assertNames(
        names(module_results),
        must.include = c(dataset_col, module_col, value_col)
      )
    } else {
      checkmate::assertList(module_results, types = "matrix")
      checkmate::assertTRUE(all(purrr::map_lgl(module_results, \(mat) {
        checkmate::checkMatrix(
          mat,
          mode = "numeric",
          row.names = "named",
          col.names = "named"
        )
      })))
    }

    # function body
    data <- if (rbh_type == "set") {
      split(
        module_results %>% dplyr::select(!!module_col, !!value_col),
        module_results[, ..dataset_col]
      ) %>%
        purrr::map(
          .,
          ~ {
            df <- .
            split(unlist(df[, ..value_col]), unlist(df[, ..module_col]))
          }
        )
    } else {
      module_results
    }

    # Finalise object
    S7::new_object(
      S7::S7_object(),
      module_data = data,
      rbh_edge_df = data.table(),
      rbh_graph = NULL,
      final_results = data.table(),
      params = list(
        rbh_type = rbh_type
      )
    )
  }
)

## methods ---------------------------------------------------------------------

### getters --------------------------------------------------------------------

#' Get the RBH results
#'
#' @description Pulls out the RBH results if you ran
#' [bixverse::generate_rbh_graph()]
#'
#' @param object The underlying class [bixverse::rbh_graph()].
#'
#' @return The data.table with the RBH result if found, otherwise NULL.
#'
#' @export
get_rbh_res <- S7::new_generic(
  name = "get_rbh_res",
  dispatch_args = "object",
  fun = function(object) {
    S7::S7_dispatch()
  }
)

#' @export
#'
#' @importFrom magrittr `%>%`
#'
#' @method get_rbh_res rbh_graph
S7::method(get_rbh_res, rbh_graph) <- function(object) {
  # Checks
  checkmate::assertClass(object, "bixverse::rbh_graph")
  # Get the data
  edge_dt <- S7::prop(object, "rbh_edge_df")
  if (is.null(edge_dt)) {
    warning("No RBH results were found. Returning NULL.")
  }
  edge_dt
}

### print ----------------------------------------------------------------------

# TODO

# snf_graph --------------------------------------------------------------------

#' Similarity network fusion
#'
#' @description
#' This class is designed to generate and store the adjacency matrices for
#' subsequent usage in similarity network fusion. There are methods to add
#' different types of adjacency matrices (pending the type of data, it can be
#' continuous, categorical and/or mixed). Subsequently, there are methods to
#' generate the final network.
#'
#' @section Properties:
#' \describe{
#'   \item{adj_matrices}{A list containing the processed adjacency matrices.}
#'   \item{snf_adj}{Matrix. The final adjacency of the fused network.}
#'   \item{params}{A (nested) list that will store all the parameters of the
#'   applied function.}
#'   \item{final_results}{data.table. Contains final results.}
#' }
#'
#' @param data Optional data to already transform into adjacency data. Any data
#' supplied will be assumed to be samples x features. The provided data taype
#' can be a data.table (for categorical and/or mixed types) or a matrix (for
#' continous types). If you provide a data.table, the function will assume the
#' first column are the sample identifiers. Please ensure that setting.
#' @param data_name Optional string. Name of the data modality.
#' @param snf_params List. The SNF parameters, see [bixverse::params_snf()]. The
#' list contains the following elements:
#' \itemize{
#'   \item k - Integer. Number of neighbours to consider.
#'   \item t - Integer. Number of iterations for the SNF algorithm.
#'   \item mu - Float. Normalisation factor for the Gaussian kernel width.
#'   \item alpha - Float. Normalisation parameter controlling the fusion
#'   strength.
#'   \item normalise - Boolean. Shall continuous values be Z-scored.
#'   \item distance_metric - String. One of
#'   `c("euclidean", "manhattan", "canberra", "cosine")`. Which distance metric
#'   to use for the continuous calculations. In case of pure categorical,
#'   Hamming will be used, for mixed data types Gower distance is used.
#' }
#' The parameters will be internally stored for subsequent usage in other
#' functions.
#'
#' @return Returns the `snf` class for further operations.
#'
#' @export
snf <- S7::new_class(
  # Names, parents
  name = "snf",
  parent = bixverse_base_class,

  # Properties, i.e., slots
  properties = list(
    adj_matrices = S7::class_list,
    snf_adj = S7::class_any,
    final_results = S7::class_data.frame,
    params = S7::class_list
  ),

  constructor = function(
    data = NULL,
    data_name = NULL,
    snf_params = params_snf()
  ) {
    # checks
    checkmate::assert(
      checkmate::testNull(data),
      checkmate::testDataTable(data),
      checkmate::testMatrix(data, mode = "numeric")
    )
    assertSNFParams(snf_params)
    checkmate::qassert(data_name, c("0", "S1"))

    if (!is.null(data)) {
      checkmate::qassert(data_name, c("S1"))
    }

    affinity <- switch(
      class(data)[1],
      "NULL" = NULL,
      "data.table" = with(
        snf_params,
        snf_process_aff_cat_mixed(data = data, k = k, mu = mu)
      ),
      "matrix" = with(
        snf_params,
        snf_process_aff_continuous(
          data = data,
          k = k,
          mu = mu,
          distance_metric = distance_metric,
          normalise = normalise
        )
      ),
      NULL
    )

    adj_matrices <- list()

    if (!is.null(affinity)) {
      adj_matrices[[data_name]] <- affinity
    }

    sample_no <- if (is.null(affinity)) NULL else nrow(affinity)

    # Finalise object
    S7::new_object(
      S7::S7_object(),
      adj_matrices = adj_matrices,
      snf_adj = NULL,
      params = list(
        "snf" = snf_params,
        "no_samples" = sample_no
      ),
      final_results = data.table()
    )
  }
)

## methods ---------------------------------------------------------------------

### getters --------------------------------------------------------------------

#' Get the SNF params
#'
#' @param object The underlying class [bixverse::snf()].
#'
#' @return Returns the stored SNF params
#'
#' @export
get_snf_params <- S7::new_generic(
  name = "get_snf_params",
  dispatch_args = "object",
  fun = function(object) {
    S7::S7_dispatch()
  }
)

#' @export
#'
#' @importFrom magrittr `%>%`
#'
#' @method get_snf_params snf
S7::method(get_snf_params, snf) <- function(object) {
  # Checks
  checkmate::assertClass(object, "bixverse::snf")
  # Get the data
  snf_params <- S7::prop(object, "params")[["snf"]]

  return(snf_params)
}


#' Get an individual affinity matrix
#'
#' @param object The underlying class [bixverse::snf()].
#' @param name String. The name of the individual data modality affinity
#' matrix to return.
#'
#' @return Returns adjcacency matrix if found.
#'
#' @export
get_snf_adjcacency_mat <- S7::new_generic(
  name = "get_snf_adjcacency_mat",
  dispatch_args = "object",
  fun = function(object, name) {
    S7::S7_dispatch()
  }
)

#' @export
#'
#' @importFrom magrittr `%>%`
#'
#' @method get_snf_adjcacency_mat snf
S7::method(get_snf_adjcacency_mat, snf) <- function(object, name) {
  # Checks
  checkmate::assertClass(object, "bixverse::snf")
  checkmate::qassert(name, "S1")
  # Get the data
  adj_matrix <- S7::prop(object, "adj_matrices")[[name]]

  if (is.null(adj_matrix)) {
    warning(paste("No adjcacency matrix with this name found. Returning NULL"))
  }

  return(adj_matrix)
}


#' Get the final SNF matrix
#'
#' @param object The underlying class [bixverse::snf()].
#'
#' @return Returns the SNF adjacency/similarity matrix.
#'
#' @export
get_snf_final_mat <- S7::new_generic(
  name = "get_snf_final_mat",
  dispatch_args = "object",
  fun = function(object) {
    S7::S7_dispatch()
  }
)

#' @export
#'
#' @importFrom magrittr `%>%`
#'
#' @method get_snf_final_mat snf
S7::method(get_snf_final_mat, snf) <- function(object) {
  # Checks
  checkmate::assertClass(object, "bixverse::snf")
  # Get the data
  res <- S7::prop(object, "snf_adj")

  if (is.null(res)) {
    warning(paste("No adjcacency matrix with this name found. Returning NULL"))
  }

  return(res)
}
