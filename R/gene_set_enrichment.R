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
  . <- `:=` <- pvals <- fdr <- hits <- target_set_lengths <- NULL
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
  . <- `:=` <- pvals <- fdr <- target_set_name <- hits <- NULL
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

#' Bixverse implementation of the simple fgsea algorithm
#'
#' @description
#' Rust-based version of the traditional permutation-based GSEA algorithm.
#'
#' @param stats Named numeric vector. The gene level statistic.
#' @param pathways List. A named list with each element containing the genes for
#' this pathway.
#' @param nperm Integer. Number of permutation tests. Defaults to `2000L`.
#' @param gsea_params List. The GSEA parameters, see [bixverse::params_gsea()]
#' wrapper function. This function generates a list containing:
#' \itemize{
#'  \item min_size - Integer. Minimum size for the gene sets.
#'  \item max_size - Integer. Maximum size for the gene sets.
#'  \item gsea_param - Float. The GSEA parameter. Defaults to `1.0`.
#' }
#' @param seed Random seed for reproducibility.
#'
#' @returns To be written.
#'
#' @export
calc_gsea_traditional = function(
  stats,
  pathways,
  nperm = 2000L,
  gsea_params = params_gsea(),
  seed = 123L
) {
  # Scope checks
  . <- `:=` <- NULL

  # Checks
  checkmate::assertNumeric(stats, min.len = 3L, finite = TRUE)
  checkmate::assertNames(names(stats))
  checkmate::assertList(pathways, types = "character")
  checkmate::assertNames(names(pathways))
  bixverse::assertGSEAParams(gsea_params)

  c(stats, pathways_clean, pathway_sizes) %<-%
    with(
      gsea_params,
      prep_stats_pathways(
        stats = stats,
        pathways = pathways,
        min_size = min_size,
        max_size = max_size
      )
    )

  gsea_stat_res <- with(
    gsea_params,
    do.call(
      rbind,
      lapply(
        pathways_clean,
        rs_calc_gsea_stats,
        stats = stats,
        gsea_param = gsea_param,
        return_leading_edge = TRUE
      )
    )
  )

  leading_edges <- mapply(
    "[",
    list(names(stats)),
    gsea_stat_res[, "leading_edge"],
    SIMPLIFY = FALSE
  )

  pathway_scores <- unlist(gsea_stat_res[, "es"])

  permutations_res_traditional <- rs_calc_gsea_stat_traditional_batch(
    stats = stats,
    pathway_scores = pathway_scores,
    pathway_sizes = as.integer(pathway_sizes),
    iters = nperm,
    seed = seed
  ) %>%
    data.table::setDT() %>%
    .[, `:=`(
      pathway_name = rownames(gsea_stat_res),
      leading_edge = leading_edges
    )]

  return(permutations_res_traditional)
}

#### simple implementation -----------------------------------------------------

#' Bixverse implementation of the simple fgsea algorithm
#'
#' @description
#' Rust-based version of the fgsea simple algorithm.
#'
#' @param stats Named numeric vector. The gene level statistic.
#' @param pathways List. A named list with each element containing the genes for
#' this pathway.
#' @param nperm Integer. Number of permutation tests. Defaults to `2000L`.
#' @param gsea_params List. The GSEA parameters, see [bixverse::params_gsea()]
#' wrapper function. This function generates a list containing:
#' \itemize{
#'  \item min_size - Integer. Minimum size for the gene sets.
#'  \item max_size - Integer. Maximum size for the gene sets.
#'  \item gsea_param - Float. The GSEA parameter. Defaults to `1.0`.
#' }
#' @param seed Random seed for reproducibility.
#'
#' @returns To be written.
#'
#' @export
calc_fgsea_simple = function(
  stats,
  pathways,
  nperm = 2000L,
  gsea_params = params_gsea(),
  seed = 123L
) {
  # Scope checks
  . <- `:=` <- NULL

  # Checks
  checkmate::assertNumeric(stats, min.len = 3L, finite = TRUE)
  checkmate::assertNames(names(stats))
  checkmate::assertList(pathways, types = "character")
  checkmate::assertNames(names(pathways))
  bixverse::assertGSEAParams(gsea_params)

  c(stats, pathways_clean, pathway_sizes) %<-%
    with(
      gsea_params,
      prep_stats_pathways(
        stats = stats,
        pathways = pathways,
        min_size = min_size,
        max_size = max_size
      )
    )

  gsea_stat_res <- with(
    gsea_params,
    do.call(
      rbind,
      lapply(
        pathways_clean,
        rs_calc_gsea_stats,
        stats = stats,
        gsea_param = gsea_param,
        return_leading_edge = TRUE
      )
    )
  )

  leading_edges <- mapply(
    "[",
    list(names(stats)),
    gsea_stat_res[, "leading_edge"],
    SIMPLIFY = FALSE
  )

  pathway_scores <- unlist(gsea_stat_res[, "es"])

  permutations_res_simple <- with(
    gsea_params,
    rs_calc_gsea_stat_cumulative_batch(
      stats = stats,
      pathway_scores = pathway_scores,
      pathway_sizes = as.integer(pathway_sizes),
      iters = nperm,
      seed = seed,
      gsea_param = gsea_param
    )
  ) %>%
    data.table::setDT() %>%
    .[, `:=`(
      pathway_name = rownames(gsea_stat_res),
      leading_edge = leading_edges
    )]

  return(permutations_res_simple)
}

#### multi level implementation ------------------------------------------------

# TODO Needs to be implemented

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
