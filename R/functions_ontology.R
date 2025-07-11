# individual functions ---------------------------------------------------------

## similarity calculations -----------------------------------------------------

#' Calculate the Resnik or Lin semantic similarity
#'
#' @description This function calculates the semantic similarities based on
#' Resnik or Lin similarity for a given ontology. Has also the option to
#' calculate a combined version of the two.
#'
#' @param terms Vector of strings. The terms in the ontology you wish to screen.
#' @param similarity_type String. One of `c("resnik", "lin", "combined")`.
#' @param ancestor_list List. Names being the term and the elements in the
#' list the names of the ancestors, see [bixverse::get_ontology_ancestry()].
#' @param ic_list List. The names being the term and the elements the
#' information content of this given term. Needs to be a single float! See
#' [bixverse::calculate_information_content()].
#'
#' @return The symmetric similarity matrix for the specified terms and semantic
#' similarity measure you chose.
#'
#' @export
#'
#' @import data.table
calculate_semantic_sim <- function(
  terms,
  similarity_type,
  ancestor_list,
  ic_list
) {
  # Checks
  checkmate::qassert(terms, "S+")
  checkmate::assertChoice(similarity_type, c("resnik", "lin", "combined"))
  checkmate::assertList(ancestor_list, types = "character")
  checkmate::assertNames(names(ancestor_list), must.include = terms)
  checkmate::assertList(ic_list, types = "double")
  checkmate::assertNames(names(ic_list), must.include = terms)

  onto_similarities <- rs_onto_semantic_sim(
    terms = terms,
    sim_type = similarity_type,
    ancestor_list = ancestor_list,
    ic_list = ic_list
  )

  # Using this one to deal with this
  matrix <- rs_upper_triangle_to_dense(
    cor_vector = onto_similarities,
    shift = 1L,
    n = length(terms)
  )
  colnames(matrix) <- rownames(matrix) <- terms

  return(matrix)
}


#' Calculate the Wang similarity matrix
#'
#' @description This function calculates the Wang similarity, based on the DAG
#' for a given ontology. The current implementation just allows for a single w
#' for the relationships.
#'
#' @param parent_child_dt data.table. The data.table with column parent and
#' child. You also need to have a type column for the Wang similarity to provide
#' the weights for the relationships.
#' @param weights Named numeric. The relationship of type to weight for this
#' specific edge. For example `c("part_of" = 0.8, "is_a" = 0.6)`.
#'
#' @return The symmetric Wang similarity matrix.
#'
#' @export
#'
#' @import data.table
calculate_wang_sim <- function(parent_child_dt, weights) {
  # Scope
  weight <- type <- NULL

  # Checks
  checkmate::assertDataTable(parent_child_dt)
  checkmate::assert(all(
    c("parent", "child", "type") %in% colnames(parent_child_dt)
  ))
  checkmate::assertNumeric(weights, min.len = 1L, names = "named")
  checkmate::assertTRUE(all(unique(parent_child_dt$type) %in% names(weights)))

  parent_child_dt <- data.table::copy(parent_child_dt)[,
    weight := weights[type]
  ]

  sim_wang <- rs_onto_sim_wang(
    parents = parent_child_dt$parent,
    children = parent_child_dt$child,
    w = parent_child_dt$weight,
    flat_matrix = FALSE
  )

  sim_wang_mat <- sim_wang$sim_mat %>%
    `colnames<-`(sim_wang$names) %>%
    `rownames<-`(sim_wang$names)

  return(sim_wang_mat)
}

## helpers ---------------------------------------------------------------------

#' Return ancestry terms from an ontology
#'
#' @description This function will return all ancestors and descendants based on
#' a provided data.table with parent-child terms
#'
#' @param parent_child_dt data.table. The data.table with column parent and
#' child.
#'
#' @return A list with
#' \itemize{
#'  \item ancestors A list with all ancestor terms.
#'  \item descendants A list with all descendant terms.
#' }
#'
#' @export
get_ontology_ancestry <- function(parent_child_dt) {
  . <- NULL

  checkmate::assertDataTable(parent_child_dt)
  checkmate::assert(all(c("parent", "child") %in% colnames(parent_child_dt)))

  # Deep copy to avoid side effects
  edge_df <- data.table::copy(parent_child_dt)
  data.table::setnames(edge_df, c("parent", "child"), c("to", "from"))

  graph <- igraph::graph_from_data_frame(edge_df[, c("from", "to")])
  ancestor_DT <- graph %>%
    igraph::ego(order = igraph::vcount(graph), mode = "out") %>%
    stats::setNames(igraph::V(graph)$name) %>%
    Map(f = names) %>%
    stack() %>%
    rev() %>%
    stats::setNames(c("from", "to")) %>%
    data.table::as.data.table() %>%
    .[, lapply(.SD, as.character)]

  ancestors <- split(ancestor_DT$to, ancestor_DT$from)
  descendants <- split(ancestor_DT$from, ancestor_DT$to)

  return(list(ancestors = ancestors, descandants = descendants))
}


#' Calculate the information content for each ontology term
#'
#' @description this function will calculate the information content of each
#' provided term based on a list of ancestors, which is a named list of terms
#' with ancestor identifiers as their values. Can be calculated using
#' [bixverse::get_ontology_ancestry()]. The information content is calculated
#' as `-log2(number descendant/total terms in the ontology)`. More information
#' can be found [here](https://yulab-smu.top/biomedical-knowledge-mining-book/semantic-similarity-overview.html).
#'
#' @param ancestor_list List. Named list of terms with ancestor identifiers as
#' their values
#'
#' @return A named list of each term and their information content as values.
#'
#' @export
#'
#' @import data.table
#' @importFrom magrittr `%>%`
calculate_information_content <- function(ancestor_list) {
  checkmate::assertList(ancestor_list)
  checkmate::assertNamed(ancestor_list)

  tab <- purrr::map_dbl(ancestor_list, length)
  information_content <- stats::setNames(
    -log(as.integer(tab) / length(ancestor_list)),
    nm = names(ancestor_list)
  )
  information_content <- as.list(information_content)

  return(information_content)
}


#' Calculates the critical value
#'
#' @description This function calculates the critical value for a given ontology
#' similarity matrix.
#'
#' @param x Numerical matrix or `ontology class`, see [bixverse::ontology()].
#' This function tends to be slower on matrices compared to `ontology class`.
#' @param alpha Float. The alpha value. For example, 0.001 would mean that the
#' critical value is smaller than 0.1 percentile of the random permutations.
#' @param permutations Number of random permutations.
#' @param seed Integer. For reproducibility purposes
#'
#' @return The critical value.
#'
#' @export
calculate_critical_value <- function(
  x,
  alpha,
  permutations = 100000L,
  seed = 10101L
) {
  # checks
  checkmate::assert(
    checkmate::test_matrix(x, mode = "numeric"),
    checkmate::test_class(x, "bixverse::ontology")
  )
  checkmate::qassert(alpha, "N1(0, 1)")
  checkmate::qassert(permutations, "I1")
  checkmate::qassert(seed, "I1")

  crit_val <- if (checkmate::test_matrix(x)) {
    rs_critval_mat(mat = x, iters = permutations, alpha = alpha, seed = seed)
  } else {
    sim_mat <- S7::prop(x, "sim_mat")
    if (is.null(sim_mat)) {
      warning("No sim_mat found in the object. Returning NULL.")
      return(NULL)
    } else {
      data <- sim_mat$get_data()[["data"]]
      rs_critval(
        values = data,
        iters = permutations,
        alpha = alpha,
        seed = seed
      )
    }
  }

  return(crit_val)
}
