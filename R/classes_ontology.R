# gene ontology gene set enrichment --------------------------------------------

#' Gene Ontology data
#'
#' @description
#' This class is used to store the gene ontology information for usage in GSE
#' elimination methods.
#'
#' @param go_data_dt A data.table that contains the gene ontology information.
#' This can be extract with for example [bixverse::get_go_human_data()].
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
#' @return Returns the class for subsequent usage.
#'
#' @export
gene_ontology_data <- S7::new_class(
  # Names, parents
  name = "gene_ontology_data",

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

#' @name print.gene_ontology_data
#'
#' @title print Method for gene_ontology_data object
#'
#' @description
#' Print a gene_ontology_data object.
#'
#' @param x An object of class `gene_ontology_data`, see
#' [bixverse::gene_ontology_data()].
#' @param ... Additional arguments (currently not used).
#'
#' @returns Invisibly returns `x`.
#'
#' @method print gene_ontology_data
S7::method(print, gene_ontology_data) <- function(x, ...) {
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


# ontology class ---------------------------------------------------------------

#' Ontology class
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
#' @return Returns the class for subsequent usage.
#'
#' @export
ontology <- S7::new_class(
  # Names, parents
  parent = bixverse_base_class,
  name = "ontology",

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

#' @name print.ontology
#'
#' @title print Method for ontology object
#'
#' @description
#' Print a ontology object.
#'
#' @param x An object of class `ontology`, see [bixverse::ontology()].
#' @param ... Additional arguments (currently not used).
#'
#' @returns Invisibly returns `x`.
#'
#' @method print ontology
S7::method(print, ontology) <- function(x, ...) {
  # Get necessary parameters
  ontology_size <- S7::prop(x, "params")[["ontology_data"]][["total_size"]]
  sim_mat <- S7::prop(x, "sim_mat")
  semantic_calculated <- ifelse(is.null(sim_mat), "No.", "Yes.")

  cat(paste(
    "Ontology class:",
    sprintf(" Size of the ontology: %i.", ontology_size),
    sprintf(" Semantic similarities calculated: %s", semantic_calculated),
    sep = "\n"
  ))

  invisible(x)
}

## getters ---------------------------------------------------------------------

#' Get the similarity matrix
#'
#' @param object `ontology class`. See [bixverse::ontology()].
#' @param as_data_table Boolean. Shall the data be returned as a long
#' data.table.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @return Returns the semantic similarity data.table from the class
#'
#' @export
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
#' @importFrom magrittr `%>%`
#'
#' @method get_sim_matrix ontology
S7::method(get_sim_matrix, ontology) <-
  function(object, as_data_table = FALSE, .verbose = TRUE) {
    checkmate::assertClass(object, "bixverse::ontology")
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
