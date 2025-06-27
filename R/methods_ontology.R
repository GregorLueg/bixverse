# methods ----------------------------------------------------------------------

## pre-processing --------------------------------------------------------------

#' Pre-process data for subsequent ontology similarity
#'
#' @description This function calculates needed information for semantic
#' similiary calculations
#'
#' @param object `ontology class`. See [bixverse::ontology()].
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @return The class with added pre-processed data for semantic similarities to
#' the properties.
#'
#' @export
pre_process_sim_onto <- S7::new_generic(
  name = "pre_process_sim_onto",
  dispatch_args = "object",
  fun = function(object, .verbose = TRUE) {
    S7::S7_dispatch()
  }
)

#' @export
#'
#' @import data.table
#' @importFrom magrittr `%>%`
#'
#' @method pre_process_sim_onto ontology
S7::method(pre_process_sim_onto, ontology) <- function(
  object,
  .verbose = TRUE
) {
  # Checks
  checkmate::assertClass(object, "bixverse::ontology")
  checkmate::qassert(.verbose, "B1")

  parent_child_dt <- S7::prop(object, "parent_child_dt")

  if (.verbose) {
    message("Identifying the ancestors in the ontology.")
  }
  c(ancestors, descendants) %<-% get_ontology_ancestry(parent_child_dt)
  if (.verbose) {
    message("Calculating the information content of each term")
  }
  information_content <- calculate_information_content(descendants)

  S7::prop(object, "outputs")[['ancestors']] <- ancestors
  S7::prop(object, "outputs")[['descendants']] <- descendants
  S7::prop(object, "outputs")[['information_content']] <- information_content

  return(object)
}

## similarities ----------------------------------------------------------------

### semantic similarities ------------------------------------------------------

#' Calculate the Resnik or Lin semantic similarity for an ontology.
#'
#' @description This function calculates the specified semantic similarities for
#' the whole ontology and adds it to the class.
#'
#' @param object `ontology class`. See [bixverse::ontology()].
#' @param sim_type String. One of `c("resnik", "lin", "combined")`.
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @return The class with added semantic similarities to the properties.
#'
#' @export
calculate_semantic_sim_onto <- S7::new_generic(
  name = "calculate_semantic_sim_onto",
  dispatch_args = "object",
  fun = function(
    object,
    sim_type = c("resnik", "lin", "combined"),
    .verbose = TRUE
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
    sim_type = c("resnik", "lin", "combined"),
    .verbose = TRUE
  ) {
    sim_type <- match.arg(sim_type)

    # Checks
    checkmate::assertClass(object, "bixverse::ontology")
    checkmate::assertChoice(sim_type, c("resnik", "lin", "combined"))
    checkmate::qassert(.verbose, "B1")

    if (is.null(S7::prop(object, "outputs")[['ancestors']])) {
      warning(
        paste(
          "No pre-processed data for semantic similarity calculation found.",
          "Running the function now."
        )
      )
      object <- pre_process_sim_onto(object = object, .verbose = .verbose)
    }

    ancestor_list <- S7::prop(object, "outputs")[['ancestors']]
    information_content_list <- S7::prop(object, "outputs")[[
      'information_content'
    ]]
    terms <- names(ancestor_list)

    if (.verbose) {
      message(sprintf(
        "Calculating the semantic similarities with type: %s.",
        sim_type
      ))
    }

    similarities <- rs_onto_semantic_sim(
      terms = terms,
      sim_type = sim_type,
      ancestor_list = ancestor_list,
      ic_list = information_content_list
    )

    final_sim <- upper_triangular_cor_mat$new(
      cor_coef = similarities,
      shift = 1L,
      features = terms
    )

    params <- list(
      sim_type = sim_type
    )

    S7::prop(object, "sim_mat") <- final_sim
    S7::prop(object, "params")[["similarity"]] <- params

    return(object)
  }

### wang similarity ------------------------------------------------------------

#' Calculate the Wang similarity for an ontology.
#'
#' @description This function calculates the Wang similarity for the whole
#' ontology and adds it to the class. The current implementation just allows for
#' a single w for the relationships.
#'
#' @param object `ontology class`. See [bixverse::ontology()].
#' @param w Float. The value assigned to the respective edge. Defaults to `0.8`.
#' @param .verbose Boolean. Controls the verbosity of the function.
#'
#' @return The class with added semantic similarities to the properties.
#'
#' @export
calculate_wang_sim_onto <- S7::new_generic(
  name = "calculate_wang_sim_onto",
  dispatch_args = "object",
  fun = function(
    object,
    w = 0.8,
    .verbose = TRUE
  ) {
    S7::S7_dispatch()
  }
)


#' @export
#'
#' @import data.table
#' @importFrom magrittr `%>%`
#' @importFrom zeallot `%<-%`
#'
#' @method calculate_wang_sim_onto ontology
S7::method(calculate_wang_sim_onto, ontology) <- function(
  object,
  w = 0.8,
  .verbose = TRUE
) {
  # checks
  checkmate::assertClass(object, "bixverse::ontology")
  checkmate::qassert(.verbose, "B1")

  # body
  parent_child_dt <- S7::prop(object, "parent_child_dt")

  c(sim_data, features) %<-%
    rs_onto_sim_wang(
      parent_child_dt$parent,
      parent_child_dt$child,
      w = w,
      flat_matrix = TRUE
    )

  final_sim <- upper_triangular_cor_mat$new(
    cor_coef = sim_data,
    shift = 1L,
    features = features
  )

  params <- list(
    sim_type = "wang",
    w = w
  )

  S7::prop(object, "sim_mat") <- final_sim
  S7::prop(object, "params")[["similarity"]] <- params

  return(object)
}

# individual functions ---------------------------------------------------------

## main functions --------------------------------------------------------------

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
  diag(matrix) <- 0
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
#' child.
#' @param w Float. The value assigned to the respective edge. Defaults to `0.8`.
#'
#' @return The symmetric Wang similarity matrix.
#'
#' @export
#'
#' @import data.table
calculate_wang_sim <- function(parent_child_dt, w = 0.8) {
  # Checks
  checkmate::assertDataTable(parent_child_dt)
  checkmate::assert(all(c("parent", "child") %in% colnames(parent_child_dt)))
  checkmate::qassert(w, "R1(0,1)")

  sim_wang <- rs_onto_sim_wang(
    parent_child_dt$parent,
    parent_child_dt$child,
    w = w,
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
#' similarity matrix
#'
#' @param x Numerical matrix or `ontology class`, see [bixverse::ontology()].
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

  data <- if (checkmate::test_matrix(x)) {
    rs_dense_to_upper_triangle(x, 1L)
  } else {
    sim_mat <- S7::prop(x, "sim_mat")
    if (is.null(sim_mat)) {
      return(NULL)
    } else {
      sim_mat$get_cor_vector()[["cor_data"]]
    }
  }

  crit_val <- rs_critval(
    values = data,
    iters = permutations,
    alpha = alpha,
    seed = seed
  )

  return(crit_val)
}
