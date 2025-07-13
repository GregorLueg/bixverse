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
  checkmate::assertList(target_genes_list, types = "character", names = "named")
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

## simplify --------------------------------------------------------------------

#' Simplify gene set results via ontologies
#'
#' @description
#' This function provides an interface to simplify overenrichment results
#' based on ontological information of the pathway origin (typical use case
#' is to simplify gene ontology results). To do so, the function will
#' calculate the Wang similarity and keep within a set of highly similar terms
#' the one with the lowest fdr. Should there be terms with the same fdr, the
#' function will keep the most specific term within the ontology.
#'
#' @param res data.table. The enrichment results. Needs to have the columns
#' `c("gene_set_name", "fdr")`.
#' @param parent_child_dt data.table. The data.table with column parent and
#' child. You also need to have a type column for the Wang similarity to provide
#' the weights for the relationships.
#' @param weights Named numeric. The relationship of type to weight for this
#' specific edge. For example `c("part_of" = 0.8, "is_a" = 0.6)`.
#' @param min_sim Float between 0 and 1. The minimum similarity that the terms
#' need to have.
#'
#' @return data.table with enrichment results.
#'
#' @export
#'
#' @importFrom magrittr `%>%`
#' @importFrom magrittr `%$%`
#' @import data.table
simplify_hypergeom_res <- function(
  res,
  parent_child_dt,
  weights,
  min_sim = 0.7
) {
  # checks
  checkmate::assertDataTable(res)
  checkmate::assertNames(names(res), must.include = c("gene_set_name", "fdr"))
  checkmate::assertDataTable(parent_child_dt)
  checkmate::assert(all(
    c("parent", "child", "type") %in% colnames(parent_child_dt)
  ))
  checkmate::assertNumeric(weights, min.len = 1L, names = "named")
  checkmate::assertTRUE(all(unique(parent_child_dt$type) %in% names(weights)))
  all_terms <- unique(c(parent_child_dt$parent, parent_child_dt$child))
  checkmate::assertTRUE(all(res$gene_set_name %in% all_terms))
  checkmate::qassert(min_sim, "N1[0, 1]")

  # function body
  ancestry <- get_ontology_ancestry(go_parent_child_dt)

  descendants <- ancestry$descandants

  wang_sims <- calculate_wang_sim(
    terms = res$gene_set_name,
    parent_child_dt = go_parent_child_dt,
    weights = weights,
    add_self = TRUE
  )

  res_combined <- merge(
    res,
    wang_sims,
    by.x = "gene_set_name",
    by.y = "term1"
  )

  to_remove <- c()

  go_ids <- unique(res_combined$gene_set_name)

  for (i in seq_along(go_ids)) {
    id <- go_ids[i]
    subset <- res_combined[term2 == id & sims >= min_sim]

    if (nrow(subset) == 1) {
      next
    }

    to_select <- which(subset$fdr == min(subset$fdr))
    if (length(to_select) == 1) {
      # case where we can just summarise by FDR
      to_remove <- append(to_remove, subset$gene_set_name[-to_select])
    } else {
      # in this case, we will keep the most specific term
      no_descendants <- purrr::map_dbl(
        descendants[subset$gene_set_name],
        length
      )
      to_select_desc <- which(no_descendants == min(no_descendants))
      to_remove <- append(to_remove, subset$gene_set_name[-to_select])
    }
  }

  res[!gene_set_name %in% to_remove]
}

## gsea ------------------------------------------------------------------------

### main functions -------------------------------------------------------------

#### original implementations --------------------------------------------------

#' Bixverse implementation of the traditional GSEA algorithm
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
#'  \item sample_size - Integer. Number of samples to iterate through for the
#'  multi-level implementation of fgsea.
#'  \item eps - Float. Boundary for calculating the p-value. Used for the multi-
#'  level implementation of fgsea.
#' }
#' @param seed Random seed for reproducibility.
#'
#' @returns A data.table with the results from the GSEA with the following
#' columns:
#' \itemize{
#'  \item es - Float. The enrichment score for this pathway.
#'  \item nes - Float. The normalised enrichment score for this pathway.
#'  \item pvals - Float. The p-value for this pathway.
#'  \item n_more_extreme - Integer. Number of permutation that had more extreme
#'  enrichment scores than the actual.
#'  \item size - Integer. The size of the pathway.
#'  \item pathway_name - Character. The name of the pathway.
#'  \item leading_edge - List of character vectors with the leading edge genes.
#'  \item fdr - Float. The adjusted pval.
#' }
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
  . <- `:=` <- pvals <- pathways_clean <- pathway_sizes <- p.adjust <- NULL

  # Checks
  checkmate::assertNumeric(stats, min.len = 3L, finite = TRUE)
  checkmate::assertNames(names(stats))
  checkmate::assertList(pathways, types = "character")
  checkmate::assertNames(names(pathways))
  assertGSEAParams(gsea_params)

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
      leading_edge = leading_edges,
      fdr = p.adjust(pvals, method = "fdr")
    )]

  return(permutations_res_traditional)
}

#### simple implementation -----------------------------------------------------

#' Bixverse implementation of the simple fgsea algorithm
#'
#' @description
#' Rust-based version of the fgsea simple algorithm. This one is permutation-
#' based and similar to the traditional implementation, but leverages some
#' clear tricks to be way faster, see Korotkevich, et al.
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
#'  \item sample_size - Integer. Number of samples to iterate through for the
#'  multi-level implementation of fgsea.
#'  \item eps - Float. Boundary for calculating the p-value. Used for the multi-
#'  level implementation of fgsea.
#' }
#' @param seed Random seed for reproducibility.
#'
#' @returns A data.table with the results from the GSEA with the following
#' columns:
#' \itemize{
#'  \item es - Float. The enrichment score for this pathway.
#'  \item nes - Float. The normalised enrichment score for this pathway.
#'  \item pvals - Float. The p-value for this pathway.
#'  \item n_more_extreme - Integer. Number of permutation that had more extreme
#'  enrichment scores than the actual.
#'  \item size - Integer. The size of the pathway.
#'  \item pathway_name - Character. The name of the pathway.
#'  \item leading_edge - List of character vectors with the leading edge genes.
#'  \item fdr - Float. The adjusted pval.
#' }
#'
#' @export
#'
#' @references Korotkevich, et al., bioRxiv
calc_fgsea_simple = function(
  stats,
  pathways,
  nperm = 2000L,
  gsea_params = params_gsea(),
  seed = 123L
) {
  # Globals scope check
  . <- `:=` <- mode_fraction <- pvals <- log2err <- pathways_clean <- NULL
  pathway_sizes <- es <- ge_zero <- n_more_extreme <- denom_prob <- p.adjust <- NULL

  # Checks
  checkmate::assertNumeric(stats, min.len = 3L, finite = TRUE)
  checkmate::assertNames(names(stats))
  checkmate::assertList(pathways, types = "character")
  checkmate::assertNames(names(pathways))
  assertGSEAParams(gsea_params)

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
      gsea_param = gsea_param,
      return_add_stats = FALSE
    )
  ) %>%
    data.table::setDT() %>%
    .[, `:=`(
      pathway_name = rownames(gsea_stat_res),
      leading_edge = leading_edges,
      fdr = p.adjust(pvals, method = "fdr")
    )]

  return(permutations_res_simple)
}

#### multi level implementation ------------------------------------------------

#' Bixverse implementation of the fgsea algorithm
#'
#' @description
#' Rust-based version of the fgsea simple and multi-level algorithm. Initially,
#' the simple method is run. For low p-values, the multi-level method is used
#' to estimate lower p-values than possible just based on the permutations, see
#' Korotkevich, et al.
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
#'  \item sample_size - Integer. Number of samples to iterate through for the
#'  multi-level implementation of fgsea.
#'  \item eps - Float. Boundary for calculating the p-value. Used for the multi-
#'  level implementation of fgsea.
#' }
#' @param seed Random seed for reproducibility.
#'
#' @returns A data.table with the results from the GSEA with the following
#' columns:
#' \itemize{
#'  \item es - Float. The enrichment score for this pathway.
#'  \item nes - Float. The normalised enrichment score for this pathway.
#'  \item pvals - Float. The p-value for this pathway.
#'  \item n_more_extreme - Integer. Number of permutation that had more extreme
#'  enrichment scores than the actual.
#'  \item size - Integer. The size of the pathway.
#'  \item pathway_name - Character. The name of the pathway.
#'  \item leading_edge - List of character vectors with the leading edge genes.
#'  \item fdr - Float. The adjusted pval.
#' }
#'
#' @export
#'
#' @references Korotkevich, et al., bioRxiv
calc_fgsea <- function(
  stats,
  pathways,
  nperm = 1000L,
  gsea_params = params_gsea(),
  seed = 123L
) {
  # Globals scope check
  . <- `:=` <- mode_fraction <- pvals <- log2err <- pathways_clean <-
    pathway_sizes <- es <- ge_zero <- le_zero <- n_more_extreme <- denom_prob <- p.adjust <- NULL

  # Checks
  checkmate::assertNumeric(stats, min.len = 3L, finite = TRUE)
  checkmate::assertNames(names(stats))
  checkmate::assertList(pathways, types = "character")
  checkmate::assertNames(names(pathways))
  assertGSEAParams(gsea_params)

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
      gsea_param = gsea_param,
      return_add_stats = TRUE
    )
  ) %>%
    data.table::setDT() %>%
    .[, `:=`(
      pathway_name = rownames(gsea_stat_res),
      leading_edge = leading_edges,
      mode_fraction = data.table::fifelse(es >= 0, ge_zero, le_zero)
    )]

  # Calculations for the multi-level version
  rs_err_res <- with(
    gsea_params,
    rs_simple_and_multi_err(
      n_more_extreme = as.integer(permutations_res_simple$n_more_extreme),
      nperm = nperm,
      sample_size = sample_size
    )
  )

  dt_simple_gsea <- permutations_res_simple[
    rs_err_res$multi_err >= rs_err_res$simple_err
  ][,
    `:=`(
      log2err = 1 /
        log(2) *
        sqrt(trigamma(n_more_extreme + 1) - trigamma(nperm + 1))
    )
  ]

  dt_multi_level <- permutations_res_simple[
    rs_err_res$multi_err < rs_err_res$simple_err
  ][, "denom_prob" := (mode_fraction + 1) / (nperm + 1)]

  rs_res = with(
    gsea_params,
    rs_calc_multi_level(
      stats = stats,
      es = dt_multi_level$es,
      pathway_size = as.integer(dt_multi_level$size),
      sample_size = sample_size,
      seed = seed,
      eps = eps,
      sign = FALSE
    )
  )

  dt_multi_level = dt_multi_level[, pvals := rs_res$pvals] %>%
    data.table::as.data.table() %>%
    .[, pvals := pmin(1, pvals / denom_prob)] %>%
    .[,
      log2err := multilevel_error(pvals, sample_size = gsea_params$sample_size)
    ] %>%
    .[, log2err := data.table::fifelse(rs_res$is_cp_ge_half, log2err, NA)]

  all_results <- list(
    dt_simple_gsea,
    dt_multi_level
  ) %>%
    data.table::rbindlist(fill = TRUE) %>%
    .[, `:=`(
      le_zero = NULL,
      ge_zero = NULL,
      mode_fraction = NULL,
      denom_prob = NULL,
      fdr = p.adjust(pvals, method = "fdr")
    )] %>%
    data.table::setorder(pvals)

  return(all_results)
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

#' Helper function to calculate the multi-level error
#'
#' @param pval Numerical vector. The p-values.
#' @param sample_size Integer. The sample size.
#'
#' @return Returns the log2error
multilevel_error <- function(pval, sample_size) {
  # Checks
  checkmate::qassert(pval, "N+")
  checkmate::qassert(sample_size, "I1")

  # Function body
  res <- sqrt(
    floor(-log2(pval) + 1) *
      (trigamma((sample_size + 1) / 2) - trigamma(sample_size + 1))
  ) /
    log(2)

  return(res)
}
