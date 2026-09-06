# gene ontology gene set enrichment --------------------------------------------

#' Gene Ontology data
#'
#' @description
#' This class is used to store the gene ontology information for usage in GSE
#' elimination methods.
#'
#' @param go_data_dt A data.table that contains the gene ontology information.
#' This can be extract with for example [bixverse::get_go_data_human()].
#' @param min_genes data.frame. Meta-data information in form of a data.frame.
#'
#' @section Properties:
#' \describe{
#'   \item{go_info}{data.table. Contains the gene ontology identifiers and
#'   names.}
#'   \item{go_to_genes}{List. Contains the genes within each gene ontology
#'   term.}
#'   \item{ancestry}{List. Contains the ancestors for each gene ontology term.}
#'   \item{levels}{List. Which gene ontology terms sit at which level.}
#'   \item{min_genes}{Integer, the minimum genes in the gene ontology term to
#'   conduct the test.}
#' }
#'
#' @returns Returns the class for subsequent usage.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # human GO restricted to terms with at least 25 genes
#' go_obj <- GeneOntologyElim(
#'   get_go_data_human(.verbose = FALSE),
#'   min_genes = 25L
#' )
#' go_obj
#' }
GeneOntologyElim <- S7::new_class(
  # Names, parents
  name = "GeneOntologyElim",

  # Properties, i.e., slots
  properties = list(
    go_info = S7::class_data.frame,
    go_to_genes = S7::class_list,
    ancestry = S7::class_list,
    levels = S7::class_list,
    min_genes = S7::class_integer
  ),
  constructor = function(go_data_dt, min_genes) {
    # Checks
    checkmate::assertDataTable(go_data_dt)
    checkmate::qassert(min_genes, "I1")
    go_data_dt <-
      copy(go_data_dt)[, `:=`(
        no_genes = purrr::map_dbl(ensembl_id, length),
        depth = sprintf("%02d", depth)
      )]

    go_data_dt <- go_data_dt[no_genes >= min_genes]

    go_info <- go_data_dt[, c("go_id", "go_name")]

    go_to_genes <- go_data_dt$ensembl_id
    names(go_to_genes) <- go_data_dt$go_id

    ancestry <- go_data_dt$ancestors
    names(ancestry) <- go_data_dt$go_id

    depth_df <- go_data_dt[, .(go_ids = list(go_id)), .(depth)]

    levels <- depth_df$go_ids
    names(levels) <- depth_df$depth

    # Finalise object
    S7::new_object(
      S7::S7_object(),
      go_info = go_info,
      go_to_genes = go_to_genes,
      ancestry = ancestry,
      levels = levels,
      min_genes = min_genes
    )
  }
)

## print ------------------------------------------------------------------------

#' @noRd
S7::method(print, GeneOntologyElim) <- function(x, ...) {
  # Get necessary parameters
  number_levels <- length(S7::prop(x, "levels"))
  number_gene_sets <- length(S7::prop(x, "go_to_genes"))
  min_genes <- S7::prop(x, "min_genes")

  cat(paste(
    "Gene ontology enrichment class:",
    sprintf(" Contains %i gene ontology terms.", number_gene_sets),
    sprintf(" Total of %i levels represented in the ontology.", number_levels),
    sprintf(" Minimum genes per term set to %i.", min_genes),
    sep = "\n"
  ))

  invisible(x)
}


# OntologySim class ------------------------------------------------------------

#' OntologySim class
#'
#' @description
#' This class is used to store any ontology and apply different methods to it.
#' Currently implemented are different types of term similarities, i.e., based
#' on semantic similarities and the Wang similarity.
#'
#' @param parent_child_dt A data.table that contains the ontological information
#' in terms of parent child relationships. Need to contain the
#' `c("parent", "child")` columns.
#' @param .verbose Boolean. Controls the verbosity of the class
#'
#' @section Properties:
#' \describe{
#'   \item{edge_dt}{data.table. Contains the parent-child relationships. (For
#'   Wang similarity also the relationship type.)}
#'   \item{outputs}{List. Contains various intermediary results used for some
#'   methods.}
#'   \item{sim_mat}{List. Contains the potentially calculated similarity
#'   matrix in form of an R6 class. Getters to access the data are provided.}
#'   \item{params}{A (nested) list that will store all the parameters of the
#'   applied function.}
#'   \item{final_results}{Final results stored in the class.}
#' }
#'
#' @returns Returns the class for subsequent usage.
#'
#' @export
#'
#' @examples
#' # wrap a toy parent-child ontology into the class
#' onto <- data.table::data.table(
#'   parent = c("a", "b", "b", "b", "c"),
#'   child = c("b", "c", "d", "e", "f"),
#'   type = c("part_of", "part_of", "part_of", "is_a", "is_a")
#' )
#' OntologySim(onto, .verbose = FALSE)
OntologySim <- S7::new_class(
  # Names, parents
  parent = BixverseBaseClass,
  name = "OntologySim",

  # Properties, i.e., slots
  properties = list(
    parent_child_dt = S7::class_data.frame,
    outputs = S7::class_list,
    sim_mat = S7::class_any,
    params = S7::class_list,
    final_results = S7::class_any
  ),
  constructor = function(parent_child_dt, .verbose = TRUE) {
    # Checks
    checkmate::assertDataTable(parent_child_dt)
    checkmate::assert(all(c("parent", "child") %in% colnames(parent_child_dt)))
    checkmate::qassert(.verbose, "B1")

    params <- list(
      ontology_data = list(
        total_size = length(unique(
          parent_child_dt$child,
          parent_child_dt$parent
        ))
      )
    )

    # Finalise object
    S7::new_object(
      S7::S7_object(),
      parent_child_dt = parent_child_dt,
      outputs = list(),
      sim_mat = NULL,
      params = params,
      final_results = NULL
    )
  }
)

## print -----------------------------------------------------------------------

#' @noRd
S7::method(print, OntologySim) <- function(x, ...) {
  # Get necessary parameters
  ontology_size <- S7::prop(x, "params")[["ontology_data"]][["total_size"]]
  sim_mat <- S7::prop(x, "sim_mat")
  semantic_calculated <- ifelse(is.null(sim_mat), "No.", "Yes.")

  cat(paste(
    "OntologySim class:",
    sprintf(" Size of the ontology: %i.", ontology_size),
    sprintf(" Semantic similarities calculated: %s", semantic_calculated),
    sep = "\n"
  ))

  invisible(x)
}

## getters ---------------------------------------------------------------------

#' Get the similarity matrix
#'
#' @param object `OntologySim class`. See [bixverse::OntologySim()].
#' @param as_data_table Boolean. Shall the data be returned as a long
#' data.table.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @returns Returns the semantic similarity data.table from the class
#'
#' @export
#'
#' @examples
#' # the Wang similarity matrix back out of the class
#' onto <- data.table::data.table(
#'   parent = c("a", "b", "b", "b", "c"),
#'   child = c("b", "c", "d", "e", "f"),
#'   type = c("part_of", "part_of", "part_of", "is_a", "is_a")
#' )
#' onto_obj <- OntologySim(onto, .verbose = FALSE)
#' onto_obj <- calculate_wang_sim_onto(
#'   onto_obj,
#'   weights = c(part_of = 0.8, is_a = 0.6),
#'   .verbose = FALSE
#' )
#' round(get_sim_matrix(onto_obj, .verbose = FALSE), 3)
get_sim_matrix <- S7::new_generic(
  name = "get_semantic_similarities",
  dispatch_args = "object",
  fun = function(object, as_data_table = FALSE, .verbose = TRUE) {
    S7::S7_dispatch()
  }
)

#' @export
#'
#' @import data.table
#' @importFrom magrittr %>%
#'
#' @method get_sim_matrix OntologySim
S7::method(get_sim_matrix, OntologySim) <-
  function(object, as_data_table = FALSE, .verbose = TRUE) {
    checkmate::assertClass(object, "bixverse::OntologySim")
    sim_mat <- S7::prop(object, "sim_mat")
    # Early return
    if (is.null(sim_mat)) {
      return(NULL)
    }
    res <- if (as_data_table) {
      sim_mat$get_data_table(.verbose = .verbose)
    } else {
      sim_mat$get_sym_matrix(.verbose = .verbose)
    }

    return(res)
  }
