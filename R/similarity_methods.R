# methods ----

# TODO Write this around an S7 class.

# helpers ----

#' Calculate the Resnik and Lin semantic similarity
#'
#' @description This function calculates the semantic similarities based on
#' Resnik and Lin similarity for a given ontology.
#'
#' @param terms Vector of strings. The terms in the ontology you wish to screen.
#' @param ancestor_list List. Names being the term and the elements in the
#' list the names of the ancestors.
#' @param ic_list List. The names being the term and the elements the
#' information content of this given term. Needs to be a single float!
#' @param max_ic Double. The maximum information content observed in the data.
#' This value will be used to calculate the normalised Resnik similarity.
#'
#' @return data.table with the following columns:
#' \itemize{
#'  \item term1 - String, the first term.
#'  \item term2 - String, the second term.
#'  \item resnik_sim - Float, the unnormalised Resnik similarity.
#'  \item lin_sim - Float, the Lin similarity.
#'  \item resnik_sim - Float, the normalised Resnik similarity (i.e., divided
#'  by max information content observed.)
#' }
#'
#' @export
#'
#' @import data.table
calculate_semantic_sim <- function(terms, ancestor_list, ic_list, max_ic) {
  # Checks
  checkmate::qassert(terms, "S+")
  checkmate::assertList(ancestor_list, types = "character")
  checkmate::assert_named(ancestor_list)
  checkmate::assertList(ic_list, types =  "double")
  checkmate::assert_named(ic_list)
  checkmate::qassert(max_ic, "N1")

  onto_similarities <- rs_onto_similarity(terms = terms,
                                          ancestor_list = ancestor_list,
                                          ic_list = ic_list) %>%
    setDT() %>%
    .[, resnik_sim_norm := resnik_sim / max_ic]

  return(bxv_sim)
}


#' Return ancestor terms from an ontology
#'
#' @description this function will return all ancestor terms based on a provided
#' data.table with parent-child terms
#'
#' @param ontology data.table. The data.table with column parent and child.
#'
#' @return A named list with ancestor terms as values
#'
#' @export
get_ontology_ancestors <- function(ontology) {
  checkmate::assertDataTable(ontology)
  checkmate::assert(all(c("parent", "child") %in% colnames(ontology)))

  # Deep copy to avoid side effects
  edge_df <- data.table::copy(ontology)
  data.table::setnames(edge_df, c("parent", "child"), c("to", "from"))

  graph = igraph::graph_from_data_frame(edge_df)
  purrr::map(igraph::V(graph), ~ {
    igraph::subcomponent(graph = graph,
                         v = .x,
                         mode = "in")$name
  })
}


#' Calculate the information content for each ontology term
#'
#' @description this function will calculate the information content of each
#' provided term based on a list of ancestors, which is a named list of terms
#' with ancestor identifiers as their values. Can be calculated using
#' [bixverse::get_ancestors()]. The information content is calculated as
#' `-log2(number descendant/total terms in the ontology)`. More information
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

  total_terms = length(unique(c(
    unlist(ancestor_list), names(ancestor_list)
  )))
  ic = as.data.table(stack(ancestor_list)) %>%
    .[, `:=`(values = as.character(values), ind = as.character(ind))] %>%
    setnames(., c("values", "ind"), c("ancestor", "descendant")) %>%
    unique() %>%
    .[, nb_desc := .N, by = ancestor] %>%
    .[, ic := -log2(nb_desc / total_terms)] %>%
    setnames(., "ancestor", "id")

  as.list(setNames(ic$ic, ic$id))
}
