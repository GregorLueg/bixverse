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
  # Scope
  ancestors <- descendants <- NULL

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

    similarities <- rs_onto_semantic_sim_mat(
      sim_type = sim_type,
      ancestor_list = ancestor_list,
      ic_list = information_content_list,
      flat_matrix = TRUE
    )

    final_sim <- upper_triangular_sym_mat$new(
      values = similarities$sim_mat,
      shift = 1L,
      features = similarities$names
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
#' ontology and adds it to the class in a memory-efficient format for
#' subsequent usage.
#'
#' @param object `ontology class`. See [bixverse::ontology()].
#' @param weights Named numeric. The relationship of type to weight for this
#' specific edge. For example `c("part_of" = 0.8, "is_a" = 0.6)`.
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
    weights,
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
  weights,
  .verbose = TRUE
) {
  # scope
  weight <- type <- NULL

  # checks
  checkmate::assertClass(object, "bixverse::ontology")
  checkmate::assertNumeric(weights, min.len = 1L, names = "named")
  checkmate::qassert(.verbose, "B1")

  # body
  parent_child_dt <- S7::prop(object, "parent_child_dt")

  # early return if no type column is found
  if (
    !checkmate::test_names(
      names(parent_child_dt),
      must.include = c("parent", "child", "type")
    )
  ) {
    warnings(paste(
      "No type column found in the column.",
      "Please consider recreating the class with the respective column.",
      "Returning class as is."
    ))
    return(object)
  }

  if (.verbose) {
    message("Calculating the Wang similarity.")
  }

  parent_child_dt <- data.table::copy(parent_child_dt)[,
    weight := weights[type]
  ]

  c(sim_data, features) %<-%
    rs_onto_sim_wang_mat(
      parents = parent_child_dt$parent,
      children = parent_child_dt$child,
      w = parent_child_dt$weight,
      flat_matrix = TRUE
    )

  final_sim <- upper_triangular_sym_mat$new(
    values = sim_data,
    shift = 1L,
    features = features
  )

  params <- list(
    sim_type = "wang",
    weights = weights
  )

  S7::prop(object, "sim_mat") <- final_sim
  S7::prop(object, "params")[["similarity"]] <- params

  return(object)
}

## filtering -------------------------------------------------------------------

#' Filter the calculated similarities
#'
#' @description This function calculates the critical value, see
#' [bixverse::calculate_critical_value()] and filters subsequently all the
#' term pairs to the ones with a value ≥ critical value.
#'
#' @param object `ontology class`. See [bixverse::ontology()].
#' @param alpha Float. The alpha value. For example, 0.001 would mean that the
#' critical value is smaller than 0.1 percentile of the random permutations.
#' @param permutations Number of random permutations.
#' @param seed Integer. For reproducibility purposes.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @return The class with filtered results added to the respective slot.
#'
#' @export
filter_similarities <- S7::new_generic(
  name = "filter_similarities",
  dispatch_args = "object",
  fun = function(
    object,
    alpha,
    permutations = 100000L,
    seed = 10101L,
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
#' @method calculate_wang_sim_onto ontology
S7::method(filter_similarities, ontology) <- function(
  object,
  alpha,
  permutations = 100000L,
  seed = 10101L,
  .verbose = TRUE
) {
  # checks
  checkmate::assertClass(object, "bixverse::ontology")
  checkmate::qassert(alpha, "N1(0, 1)")
  checkmate::qassert(permutations, "I1")
  checkmate::qassert(seed, "I1")
  checkmate::qassert(.verbose, "B1")

  # early return
  sim_mat <- S7::prop(object, "sim_mat")

  if (is.null(sim_mat)) {
    warning(paste(
      "No similarity matrix found in the class. You need to run one of the methods",
      "Returning class as is"
    ))
    return(object)
  }

  if (.verbose) {
    message(sprintf("Calculating critical value with alpha: %.4f", alpha))
  }

  critval <- calculate_critical_value(
    x = object,
    alpha = alpha,
    permutations = permutations,
    seed = seed
  )

  data <- sim_mat$get_data()

  if (.verbose) {
    message(sprintf("Filtering data with critical value %.4f", critval))
  }

  filtered_data <- data.table::setDT(rs_filter_onto_sim(
    sim_vals = data$data,
    names = data$features,
    threshold = critval
  ))

  params <- list(
    alpha = alpha,
    permutations = permutations,
    seed = seed,
    crit_val = critval
  )

  S7::prop(object, "params")[["sim_filters"]] <- params
  S7::prop(object, "final_results") <- filtered_data

  return(object)
}
