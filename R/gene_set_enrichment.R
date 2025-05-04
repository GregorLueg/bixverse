# gse functions ----------------------------------------------------------------

## hypergeom functions ---------------------------------------------------------

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

  target_genes_length <- length(target_genes)

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
      target_set_lengths := target_genes_length
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

## gsea ------------------------------------------------------------------------

### main functions -------------------------------------------------------------

#### original implementations --------------------------------------------------

#### simple implementation -----------------------------------------------------

#' Bixverse implementation of the simple fgsea algorithm
#'
#' @param stats ...
#' @param pathways ...
#' @param nperm ...
#' @param gsea_param ...
#'
#' @returns ...
#'
#' @export
fgsea_simple = function(
  stats,
  pathways,
  nperm = 2000L,
  gsea_params = params_gsea(),
  seed = 123L
) {
  # Checks
  checkmate::assertNumeric(stats, min.len = 3L, finite = TRUE)
  checkmate::assertNames(names(stats))
  checkmate::assertList(pathways, types = "character")
  checkmate::assertNames(names(pathways))
  assertGSEAParams(gsea_params)

  # Function body
  min_size <- gsea_params$min_size
  max_size <- gsea_params$max_size
  gsea_param <- gsea_params$gsea_param

  ## Prepare the data
  c(stats, pathways_clean, pathway_sizes) %<-%
    prep_stats_pathways(
      stats = stats,
      pathways = pathways,
      min_size = min_size,
      max_size = max_size
    )

  ## Calculate the actual enrichment
  gsea_stat_res <- do.call(
    rbind,
    lapply(
      pathways_clean,
      rs_calc_gsea_stats,
      stats = stats,
      gsea_param = 1,
      return_leading_edge = TRUE
    )
  )

  leading_edges <- mapply(
    "[",
    list(names(stats)),
    gsea_stat_res[, "leading_edge"],
    SIMPLIFY = FALSE
  )

  pathway_scores <- unlist(gsea_stat_res[, "gene_stat"])

  results_dt <- rs_calc_gsea_stat_cumulative_batch(
    stats = stats,
    pathway_scores = pathway_scores,
    pathway_sizes = as.integer(pathway_sizes),
    iters = nperm,
    seed = seed,
    gsea_param = 1
  ) %>%
    data.table::as.data.table() %>%
    .[, pathway := seq_along(pathway_scores)] %>%
    .[, `:=`(
      le_zero_mean = le_zero_sum / le_zero,
      ge_zero_mean = ge_zero_sum / ge_zero
    )]

  results_dt[, ES := pathway_scores[pathway]]

  results_dt[, NES := as.numeric(NA)]

  results_dt[
    (ES > 0 & ge_zero_mean != 0) | (ES <= 0 & le_zero_mean != 0),
    NES := ES / ifelse(ES > 0, ge_zero_mean, abs(le_zero_mean))
  ]

  results_dt[, pval := as.numeric(NA)]
  results_dt[
    !is.na(NES),
    pval := pmin((1 + le_es) / (1 + le_zero), (1 + ge_es) / (1 + ge_zero))
  ]

  results_dt[,
    `:=`(
      pathway = names(pathway),
      pathway_size = as.integer(pathway_sizes),
      leading_edges = leading_edges
    )
  ]

  return(results_dt)
}

#### helpers -------------------------------------------------------------------

#' Helper function to prepare data for GSEA
#'
#' @param stats Named numeric vector. The gene level statistic.
#' @param pathways List. A named list with each element containing the genes for
#' this pathway.
#' @param min_size Integer. Minimum size of the gene sets. Defaults to `5L`.
#' @param max_size Integer. Maximum size of the gene sets. Defaults to `500L`.
#'
#' @returns A list with the following elements:
#' \itemize{
#'  \item stats - The sorted gene level statistics
#'  \item pathways - The prepared pathways for further usage in the functions
#'  \item pathway_sizes - The sizes of the pathways
#' }
prep_stats_pathways <- function(
  stats,
  pathways,
  min_size = 5L,
  max_size = 500L
) {
  # Checks
  checkmate::assertNumeric(stats, min.len = 3L, finite = TRUE)
  checkmate::assertNames(names(stats))
  checkmate::qassert(min_size, "I1[3,)")
  checkmate::qassert(max_size, "I1[4,)")
  checkmate::assertList(pathways, types = "character")
  checkmate::assertNames(names(pathways))

  stats <- sort(stats, decreasing = TRUE)
  gene_universe <- names(stats)
  pathways <- rs_get_gs_indices(
    gene_universe = gene_universe,
    pathway_list = pathways
  )
  pathway_sizes <- purrr::map_dbl(pathways, length)
  to_keep <- pathway_sizes >= min_size & pathway_sizes <= max_size
  pathways <- pathways[to_keep]
  pathway_sizes <- pathway_sizes[to_keep]

  return(
    list(
      stats = stats,
      pathways = pathways,
      pathway_sizes = pathway_sizes
    )
  )
}
