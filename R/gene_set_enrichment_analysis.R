# main functions ---------------------------------------------------------------

## simple implementation -------------------------------------------------------

fgsea_simple = function(
  stats,
  pathways,
  min_size = 5L,
  max_size = 500L,
  nperm = 2000L,
  gsea_param = 1
) {
  # Checks
  checkmate::assertNumeric(stats, min.len = 3L, finite = TRUE)
  checkmate::assertNames(names(stats))
  checkmate::qassert(min_size, "I1[3,)")
  checkmate::qassert(max_size, "I1[4,)")
  checkmate::assertList(pathways, types = "character")
  checkmate::assertNames(names(pathways))
  checkmate::qassert(nperm, "I1")
  checkmate::qassert(gsea_param, "N1")

  # Function body

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
    seed = 123,
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

## helpers ---------------------------------------------------------------------

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
