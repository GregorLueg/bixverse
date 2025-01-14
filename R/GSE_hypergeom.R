# Gene Set Enrichment functions ----

#' Gene set enrichment (GSE) based on a hypergeometric test.
#'
#' @description
#' Takes a set of target genes, a list of gene sets and calculates a p-value (hypergeometric test) and odds ratio (OR) against
#' all the gene sets. Also applies a multiple hypothesis correction (BH) to the p-values.
#'
#' @param target_genes Character vector. GeneID(s) of the target genes.
#' @param gene_set_list Named list of character vectors. Names should represent the gene sets, pathways, and the elements the
#' genes within the respective gene set.
#' @param gene_universe Optional character vector. If you would like to specify specifically the gene universe. If set to NULL,
#' the function will default to all represented genes in the `gene_set_list`.
#' @param threshold Float between 0 and 1 to filter on the FDR. Default: 0.05. If NULL everything is returned.
#' @param verbose Boolean. Controls verbosity of the function.
#'
#' @return data.table with enrichment results.
#'
#' @export
#'
#' @importFrom magrittr `%>%`
#' @importFrom magrittr `%$%`
#' @import data.table
GSE_hypergeometric <- function(target_genes,
                               gene_set_list,
                               gene_universe = NULL,
                               threshold = 0.05,
                               minimum_overlap = 3L,
                               verbose = T) {
  # Input checks
  checkmate::qassert(target_genes, "S+")
  checkmate::assertList(gene_set_list, types = "character")
  checkmate::qassert(gene_universe, c("0", "S+"))
  checkmate::qassert(threshold, c("R1[0,1]", "0"))
  checkmate::qassert(minimum_overlap, "I1")
  checkmate::qassert(verbose, "B1")
  # Function body
  if (is.null(gene_universe)) {
    if (verbose) message("No gene universe given. Function will use the represented genes in the pathways/gene sets as reference.")
    gene_universe <- unique(unlist(gene_set_list))
  }

  gse_results
}
