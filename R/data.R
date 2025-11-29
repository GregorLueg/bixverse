#' Cell cycle genes
#'
#' (Human) cell cycle genes for scoring in single cell, please see
#' [bixverse::module_scores_sc()] and Tirosh et al, Science (2016). These
#' genes are based on the updated version in Seurat.
#'
#' @format ## `cell_cycle_genes`
#' A data.table
#' \describe{
#'   \item{hgnc_symbol}{The (HGNC) gene symbol}
#'   \item{ensembl_gene_id}{The corresponding ensembl identifiers.}
#'   \item{set}{To which cell cycle phase the gene belongs}
#' }
"cell_cycle_genes"
