
#' Calculate the Resnik and Lin semantic similarity
#'
#' @description This function calculates the Resnik or Lin similarity for a given ontology.
#'
#' @param terms Vector of strings. The terms in the ontology you wish to screen.
#' @param ancestor_list R list with names being the term and the elements in the list the names
#' of the ancestors.
#' @param ic_list R list with the names being the term and the elements the information content
#' of this given term. Needs to be a single float!
#' @param similarity_type String. Can take values of `c("resnik", "lin")`.
#' @param max_ic Double. The maximum information content observed in the data. This will return
#' normalised Resnik (i.e., scaled between 0 and 1). If set to 1, it will return the Resnik distance. Default = 1
#' @param return_matrix Boolean. Return the similarity as a matrix or vector, Default = TRUE
#'
#' @export
get_similarity_matrix <- function(
    terms,
    ancestor_list,
    ic_list,
    similarity_type = c("resnik", "lin"),
    return_matrix = TRUE,
    max_ic = 1){
  checkmate::assert_names(similarity_type, subset.of = c("resnik", "lin"))
  checkmate::assert_named(ancestor_list)
  checkmate::assert(all(terms %in% names(ancestor_list)))

  if(length(similarity_type) == 2){
    bxv_sim <- rs_onto_similarity_both(
      terms = terms,
      ancestor_list = ancestor_list,
      ic_list = ic_list,
      max_ic = max_ic)

    if(return_matrix){
      bxv_sim_resnik <- bixverse:::upper_triangular_cor_mat$new(
        cor_coef = bxv_sim$resnik_sim,
        features = bxv_sim$terms,
        shift = 1L
      )$get_cor_matrix()
      bxv_sim_lin <- bixverse:::upper_triangular_cor_mat$new(
        cor_coef = bxv_sim$lin_sim,
        features = bxv_sim$terms,
        shift = 1L
      )$get_cor_matrix()
      bxv_sim <- list(terms = bxv_sim$terms,
                      resnik_sim = bxv_sim_resnik,
                      lin_sim = bxv_sim_lin)
    }
  }else{
    bxv_sim <- rs_onto_similarity(
      terms = terms,
      ancestor_list = ancestor_list,
      ic_list = ic_list,
      max_ic = max_ic,
      similarity_type = similarity_type
    )

    if(return_matrix){
      bxv_sim_matrix <- bixverse:::upper_triangular_cor_mat$new(
        cor_coef = bxv_sim$similarities,
        features = bxv_sim$terms,
        shift = 1L
      )$get_cor_matrix()
      name = paste(similarity_type, "sim", sep = "_")
      bxv_sim <- list(bxv_sim$terms,
                      bxv_sim_matrix)
      names(bxv_sim) <- c("terms", name)
    }
  }
  return(bxv_sim)
}



#' Return ancestor terms from an ontology
#'
#' @description this function will return all ancestor terms based on a provided dataframe with parent-child terms
#'
#' @param ontology Dataframe. The dataframe with column parent and child.
#'
#' @return a named list with ancestor terms as values
#'
#' @export
#'
#' @import igraph
get_ancestors <- function(ontology){
  checkmate::assert(all(c("parent", "child") %in% colnames(ontology)))

  graph = igraph::graph_from_data_frame(parent_child %>% setnames(c("parent", "child"), c("to", "from")))
  purrr::map(igraph::V(graph),~{
    igraph::subcomponent(graph = graph, v = .x, mode = "in")$name
  })
}



#' Calculate the information content for each ontology term
#'
#' @description this function will calculate the information content of each provided term based on a
#' list of ancestors, which is a named list of terms with ancestor identifiers as their values. Can be
#' calculated using [bixverse::get_ancestors()].
#' The information content is calculated as -log2(number descendant/total terms in the ontology). More information
#' can be found [here](https://yulab-smu.top/biomedical-knowledge-mining-book/semantic-similarity-overview.html).
#'
#' @param ancestor_list List. named list of terms with ancestor identifiers as their values
#'
#' @return a named list of each term and their information content as values.
#'
#' @export
#'
#' @import polars
#' @import data.table
calculate_information_content <- function(ancestor_list){
  checkmate::assertList(ancestor_list)

  total_terms = length(unique(c(unlist(ancestor_list), names(ancestor_list))))
  ic = as.data.table(stack(ancestor_list)) %>%
                      .[, `:=`(values = as.character(values), ind = as.character(ind))] %>%
                      setnames(., c("values", "ind"), c("ancestor", "descendant")) %>%
    unique() %>%
    .[, nb_desc := .N, by = ancestor] %>%
    .[, ic := -log2(nb_desc/total_terms)] %>%
    setnames(., "ancestor", "id")
  as.list(setNames(ic$ic, ic$id))
}
