# methods ----

#' Calculate the Resnik and Lin semantic similarity for an ontology.
#'
#' @description This function calculates the specified semantic similarities for
#' the whole ontology, calculate a critical values based on `random_sample_no`
#' (based on a user-specified alpha) and return all the term-term pairs above
#' that threshold.
#'
#' @param object `ontology class`. See [bixverse::ontology()].
#' @param sim_type String. One of `c("resnik", "lin", "combined")`.
#' @param alpha Float. Between 0 and 1. A value of 0.01 indicates that the
#' critical value is at the 1% quantile of the random sample of similarities.
#' @param random_sample_no Integer. Number of random samples to use to estimate
#' the critical value.
#' @param seed Integer. Random seed for sampling reproducibility.
#'
#' @return The class with added semantic similarities to the properties.
#'
#' @export
calculate_semantic_sim_onto <- S7::new_generic(
  name = "calculate_semantic_sim_onto",
  dispatch_args = "object",
  fun = function(
    object,
    sim_type,
    alpha = 0.01,
    random_sample_no = 100000L,
    seed = 42L
  ) {
    S7::S7_dispatch()
  }
)

#' @export
#'
#' @import data.table
#' @importFrom magrittr `%>%`
#'
#' @method calculate_semantic_sim_onto ontology
S7::method(calculate_semantic_sim_onto, ontology) <-
  function(
    object,
    sim_type,
    alpha = 0.01,
    random_sample_no = 100000L,
    seed = 42L
  ) {
    # Checks
    checkmate::assertClass(object, "bixverse::ontology")
    checkmate::assertChoice(sim_type, c("resnik", "lin", "combined"))
    checkmate::qassert(alpha, "R1[0,1]")
    checkmate::qassert(random_sample_no, "I1")
    checkmate::qassert(seed, "I1")

    ancestor_list <- S7::prop(object, "ancestor_list")
    information_content_list <- S7::prop(object, "information_content_list")
    terms <- names(ancestor_list)

    similarities <- rs_onto_similarity_filtered(
      terms = terms,
      sim_type = sim_type,
      alpha = alpha,
      ancestor_list = ancestor_list,
      ic_list = information_content_list,
      iters = random_sample_no,
      seed = seed
    )

    params <- list(
      critval = similarities$critval,
      random_sample_no = random_sample_no,
      alpha = alpha,
      seed = seed,
      sim_type = sim_type
    )

    similarities_dt <- setDT(similarities[c("term1", "term2", "filtered_sim")])

    S7::prop(object, "semantic_similarities") <- similarities_dt
    S7::prop(object, "params")[["semantic_similarity"]] <- params

    return(object)
  }

# individual functions ---------------------------------------------------------

## main functions --------------------------------------------------------------

#' Calculate the Resnik and Lin semantic similarity
#'
#' @description This function calculates the semantic similarities based on
#' Resnik and Lin similarity for a given ontology.
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
#' similarity measure you chose. Self similarity is set to 0.
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

  onto_similarities <- rs_onto_similarity(
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
  diag(matrix) <- 0
  colnames(matrix) <- rownames(matrix) <- terms

  return(matrix)
}

## helpers ---------------------------------------------------------------------

#' Return ancestor terms from an ontology
#'
#' @description this function will return all ancestor terms based on a provided
#' data.table with parent-child terms
#'
#' @param parent_child_dt data.table. The data.table with column parent and
#' child.
#'
#' @return A named list with ancestor terms as values
#'
#' @export
get_ontology_ancestry <- function(parent_child_dt) {
  checkmate::assertDataTable(parent_child_dt)
  checkmate::assert(all(c("parent", "child") %in% colnames(parent_child_dt)))

  # Deep copy to avoid side effects
  edge_df <- data.table::copy(parent_child_dt)
  data.table::setnames(edge_df, c("parent", "child"), c("to", "from"))

  graph <- igraph::graph_from_data_frame(edge_df[, c("from", "to")])
  ancestor_DT <- graph %>%
    igraph::ego(order = igraph::vcount(graph), mode = "out") %>%
    setNames(igraph::V(graph)$name) %>%
    Map(f = names) %>%
    stack() %>%
    rev() %>%
    setNames(c("from", "to")) %>%
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
  information_content <- setNames(
    -log(as.integer(tab) / length(ancestor_list)),
    nm = names(ancestor_list)
  )
  information_content <- as.list(information_content)

  return(information_content)
}
