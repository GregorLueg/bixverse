# Gene Set Enrichment functions ----

## Hypergeometric functions ----

#' Gene set enrichment (GSE) based on a hypergeometric test.
#'
#' @description
#' Takes a set of target genes, a list of gene sets and calculates a p-value
#' (hypergeometric test) and odds ratio (OR) against all the gene sets. Also
#' applies a multiple hypothesis correction (BH) to the p-values.
#'
#' @param target_genes Character vector. GeneID(s) of the target genes.
#' @param gene_set_list Named list of character vectors. Names should represent
#' the gene sets, pathways, and the elements the genes within the respective
#' gene set.
#' @param gene_universe Optional character vector. If you would like to specify
#' specifically the gene universe. If set to NULL, the function will default to
#' all represented genes in the `gene_set_list`.
#' @param threshold Float between 0 and 1 to filter on the fdr. Default: 0.05.
#' If 1 everything is returned.
#' @param minimum_overlap Number of minimum overlap between the target genes
#' and the respective gene set.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @return data.table with enrichment results.
#'
#' @export
#'
#' @importFrom magrittr `%>%`
#' @importFrom magrittr `%$%`
#' @import data.table
gse_hypergeometric <- function(
  target_genes,
  gene_set_list,
  gene_universe = NULL,
  threshold = 0.05,
  minimum_overlap = 3L,
  .verbose = FALSE
) {
  # Avoid check issues
  `.` <- pvals <- fdr <- hits <- NULL
  # Input checks
  checkmate::qassert(target_genes, "S+")
  checkmate::assertList(gene_set_list, types = "character")
  checkmate::qassert(names(gene_set_list), "S+")
  checkmate::qassert(gene_universe, c("0", "S+"))
  checkmate::qassert(threshold, c("R1[0,1]", "0"))
  checkmate::qassert(minimum_overlap, "I1")
  checkmate::qassert(.verbose, "B1")
  # Function body
  if (is.null(gene_universe)) {
    if (.verbose) {
      message(
        "No gene universe given. Function will use the represented genes in the
        pathways/gene sets as reference."
      )
    }
    gene_universe <- unique(unlist(gene_set_list))
  }

  target_set_lengths = sapply(target_genes_list, length)

  gse_results <- rs_hypergeom_test(
    target_genes = target_genes,
    gene_sets = gene_set_list,
    gene_universe = gene_universe
  )

  gse_results <-
    data.table::data.table(do.call(cbind, gse_results)) %>%
    .[, `:=`(
      gene_set_name = names(gene_set_list),
      fdr = p.adjust(pvals, method = "BH")
    )] %>%
    data.table::setcolorder(
      .,
      c(
        "gene_set_name",
        "odds_ratios",
        "pvals",
        "fdr",
        "hits",
        "gene_set_lengths"
      )
    ) %>%
    .[(fdr <= threshold) & (hits >= minimum_overlap)] %>%
    data.table::setorder(., pvals) %>%
    .[,
      target_set_lengths := target_set_lengths[match(
        target_set_name,
        names(target_set_lengths)
      )]
    ]

  gse_results
}


#' Gene set enrichment (GSE) based on a hypergeometric test over a list.
#'
#' @description
#' Takes a set of list of target genes, a list of gene sets and calculates a
#' p-value (hypergeometric test) and odds ratio (OR) against all the gene sets.
#' Also applies a multiple hypothesis correction (BH) to the p-values.
#'
#' @param target_genes_list  Named list of character vectors. Names should
#' represent the identifiers of the target genes and the elements the genes.
#' @param gene_set_list Named list of character vectors. Names should represent
#' the gene sets, pathways, and the elements the genes within the respective
#' gene set.
#' @param gene_universe Optional character vector. If you would like to specify
#' specifically the gene universe. If set to NULL, the function will default to
#' all represented genes in the `gene_set_list`.
#' @param threshold Float between 0 and 1 to filter on the fdr. Default: 0.05.
#' If NULL everything is returned.
#' @param minimum_overlap Number of minimum overlap between the target genes and
#' the respective gene set.
#' @param .verbose Boolean. Controls verbosity of the function.
#'
#' @return data.table with enrichment results.
#'
#' @export
#'
#' @importFrom magrittr `%>%`
#' @importFrom magrittr `%$%`
#' @import data.table
gse_hypergeometric_list <- function(
  target_genes_list,
  gene_set_list,
  gene_universe = NULL,
  threshold = 0.05,
  minimum_overlap = 3L,
  .verbose = FALSE
) {
  # Avoid check issues
  `.` <- pvals <- fdr <- target_set_name <- NULL
  # Input checks
  checkmate::assertList(target_genes_list, types = "character")
  checkmate::qassert(names(target_genes_list), "S+")
  checkmate::assertList(gene_set_list, types = "character")
  checkmate::qassert(names(gene_set_list), "S+")
  checkmate::qassert(gene_universe, c("0", "S+"))
  checkmate::qassert(threshold, c("R1[0,1]", "0"))
  checkmate::qassert(minimum_overlap, "I1")
  checkmate::qassert(.verbose, "B1")
  # Function body
  if (is.null(gene_universe)) {
    if (.verbose) {
      message(
        "No gene universe given. Function will use the represented genes in the
        pathways/gene sets as reference."
      )
    }
    gene_universe <- unique(unlist(gene_set_list))
  }

  target_set_lengths = sapply(target_genes_list, length)

  gse_results <- rs_hypergeom_test_list(
    target_genes_list = target_genes_list,
    gene_sets = gene_set_list,
    gene_universe = gene_universe
  )

  gse_results <-
    data.table::data.table(do.call(cbind, gse_results)) %>%
    .[, `:=`(
      gene_set_name = rep(names(gene_set_list), length(target_genes_list)),
      target_set_name = rep(
        names(target_genes_list),
        each = length(gene_set_list)
      )
    )] %>%
    .[, fdr := p.adjust(pvals, method = "BH"), by = target_set_name] %>%
    data.table::setcolorder(
      .,
      c(
        "target_set_name",
        "odds_ratios",
        "pvals",
        "fdr",
        "hits",
        "gene_set_lengths"
      )
    ) %>%
    .[(fdr <= threshold) & (hits >= minimum_overlap)] %>%
    data.table::setorder(., pvals) %>%
    .[,
      target_set_lengths := target_set_lengths[match(
        target_set_name,
        names(target_set_lengths)
      )]
    ]

  gse_results
}
