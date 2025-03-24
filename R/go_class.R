# class ------------------------------------------------------------------------

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
      go_data_dt[, `:=`(
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

# print ------------------------------------------------------------------------

#' @name print.gene_ontology_data
#' @title print Method for gene_ontology_data object
#'
#' @description
#' Print a gene_ontology_data object.
#'
#' @param x An object of class `gene_ontology_data`.
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
    "Gene ontology class:",
    sprintf(" Contains %i gene ontology terms.", number_gene_sets),
    sprintf(" Total of %i levels represented in the ontology.", number_levels),
    sprintf(" Minimum genes per term set to %i.", min_genes),
    sep = "\n"
  ))

  invisible(x)
}

