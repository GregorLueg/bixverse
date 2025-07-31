# hypergeometric versions ------------------------------------------------------

#' Run gene ontology enrichment with elimination method.
#'
#' @description
#' This method takes the gene_ontology_data and a target gene set and performs
#' an GSE enrichment leveraging ontological information. It starts at the lowest
#' levels of the ontology and tests if there is significant enrichment for any
#' GO terms. If the threshold of the p-value is below the elimination threshold,
#' the genes from this term will be removed from all its ancestors. The function
#' then proceeds to the next level of the ontology and repeats the process. The
#' gene universe will be automatically set to every gene represented in the
#' ontology.
#'
#' @param object The underlying class, see [bixverse::gene_ontology_data()].
#' @param target_genes String. The target genes you wish to apply the GSEA over.
#' @param minimum_overlap Integer. Threshold for the minimal overlap.
#' @param fdr_threshold Float. Threshold for maximum fdr to include in the
#' output.
#' @param elim_threshold Float. Threshold from which p-value onwards the
#' elimination on the ancestors shall be conducted.
#' @param min_genes Integer. Minimum number of genes that have to be included in
#' the gene ontology term. If NULL, it will default to the number of minimum
#' genes stored in `gene_ontology_data`.
#'
#' @return data.table with enrichment results.
#'
#' @export
gse_go_elim_method <- S7::new_generic(
  name = "gse_go_elim_method",
  dispatch_args = "object",
  fun = function(
    object,
    target_genes,
    minimum_overlap = 3L,
    fdr_threshold = 0.05,
    elim_threshold = 0.05,
    min_genes = NULL
  ) {
    S7::S7_dispatch()
  }
)


#' @export
#'
#' @import data.table
#' @importFrom magrittr `%>%`
#'
#' @method gse_go_elim_method gene_ontology_data
S7::method(gse_go_elim_method, gene_ontology_data) <-
  function(
    object,
    target_genes,
    minimum_overlap = 3L,
    fdr_threshold = 0.05,
    elim_threshold = 0.05,
    min_genes = NULL
  ) {
    # Scope checks
    . <- pvals <- fdr <- hits <- NULL
    # First check
    checkmate::assertClass(object, "bixverse::gene_ontology_data")
    checkmate::qassert(target_genes, "S+")
    checkmate::qassert(fdr_threshold, "R+[0,1]")
    checkmate::qassert(elim_threshold, "R+[0,1]")
    checkmate::qassert(minimum_overlap, "I1")
    checkmate::qassert(min_genes, c("0", "I1"))
    # Extract relevant data from the S7 object
    if (is.null(min_genes)) {
      min_genes <- S7::prop(object, "min_genes")
    }
    levels <- names(S7::prop(object, "levels"))
    go_info <- S7::prop(object, "go_info")

    gene_universe_length <- length(unique(unlist(S7::prop(
      object,
      "go_to_genes"
    ))))

    rs_results_go <- rs_gse_geom_elim(
      target_genes = target_genes,
      levels = levels,
      go_obj = object,
      gene_universe_length = gene_universe_length,
      min_genes = min_genes,
      elim_threshold = elim_threshold,
      min_overlap = minimum_overlap,
      fdr_threshold = fdr_threshold
    )

    results_go_dt <- data.table(do.call(cbind, rs_results_go[-1])) %>%
      .[, go_id := rs_results_go$go_ids, ] %>%
      data.table::setcolorder(
        .,
        c(
          "go_id",
          "odds_ratios",
          "pvals",
          "fdr",
          "hits",
          "gene_set_lengths"
        )
      )

    results_go_dt <- merge(go_info, results_go_dt, by = "go_id") %>%
      data.table::setorder(., pvals)

    results_go_dt
  }


#' Run gene ontology enrichment with elimination method over a list.
#'
#' @description
#' This method takes the gene_ontology_data and a list of target gene sets and
#' performs an GSE enrichment leveraging ontological information. It starts at
#' the lowest levels of the ontology and tests if there is significant
#' enrichment for any GO terms. If the threshold of the p-value is below the
#' elimination threshold, the genes from this term will be removed from all its
#' ancestors. The function then proceeds to the next level of the ontology and
#' repeats the process. The class will leverage Rust threading to parallelise
#' the process. The gene universe will be automatically set to every gene
#' represented in the ontology.
#'
#' @param object The underlying class, see [bixverse::gene_ontology_data()].
#' @param target_gene_list List. The target genes list you wish to apply the
#' gene set enrichment analysis over.
#' @param minimum_overlap Integer. Threshold for the minimal overlap.
#' @param fdr_threshold Float. Threshold for maximum fdr to include in the
#' output.
#' @param elim_threshold Float. Threshold from which p-value onwards the
#' elimination on the ancestors shall be conducted.
#' @param min_genes Integer. Minimum number of genes that have to be included in
#' the gene ontology term. If NULL, it will default to the number of minimum
#' genes stored in `gene_ontology_data`.
#'
#' @return data.table with enrichment results.
#'
#' @export
gse_go_elim_method_list <- S7::new_generic(
  name = "gse_go_elim_method_list",
  dispatch_args = "object",
  fun = function(
    object,
    target_gene_list,
    minimum_overlap = 3L,
    fdr_threshold = 0.05,
    elim_threshold = 0.05,
    min_genes = NULL
  ) {
    S7::S7_dispatch()
  }
)


#' @export
#'
#' @import data.table
#' @importFrom magrittr `%>%`
#'
#' @method gse_go_elim_method_list gene_ontology_data
S7::method(gse_go_elim_method_list, gene_ontology_data) <-
  function(
    object,
    target_gene_list,
    minimum_overlap = 3L,
    fdr_threshold = 0.05,
    elim_threshold = 0.05,
    min_genes = NULL
  ) {
    # Scope checks
    . <- pvals <- fdr <- hits <- target_set_name <- NULL
    # First check
    checkmate::assertClass(object, "bixverse::gene_ontology_data")
    checkmate::assertList(
      target_gene_list,
      types = "character",
      names = "named"
    )
    checkmate::qassert(fdr_threshold, "R+[0,1]")
    checkmate::qassert(elim_threshold, "R+[0,1]")
    checkmate::qassert(minimum_overlap, "I1")
    checkmate::qassert(min_genes, c("0", "I1"))

    # Extract relevant data from the S7 object
    if (is.null(min_genes)) {
      min_genes <- S7::prop(object, "min_genes")
    }
    levels <- names(S7::prop(object, "levels"))
    go_info <- S7::prop(object, "go_info")

    gene_universe_length <- length(unique(unlist(S7::prop(
      object,
      "go_to_genes"
    ))))

    rs_results_go <- rs_gse_geom_elim_list(
      target_genes_list = target_gene_list,
      levels = levels,
      go_obj = object,
      gene_universe_length = gene_universe_length,
      min_genes = min_genes,
      elim_threshold = elim_threshold,
      min_overlap = minimum_overlap,
      fdr_threshold = fdr_threshold
    )

    cols_to_select <- c(
      "pvals",
      "fdr",
      "odds_ratios",
      "hits",
      "gene_set_lengths"
    )

    results_go_dt <- data.table(do.call(
      cbind,
      rs_results_go[cols_to_select]
    )) %>%
      .[, `:=`(
        go_id = rs_results_go$go_ids,
        target_set_name = rep(names(target_gene_list), rs_results_go$no_test)
      )]

    results_go_dt <- merge(results_go_dt, go_info, by = "go_id") %>%
      data.table::setorder(., pvals) %>%
      data.table::setcolorder(
        .,
        c(
          "target_set_name",
          "go_name",
          "go_id",
          "odds_ratios",
          "pvals",
          "fdr",
          "hits",
          "gene_set_lengths"
        )
      )

    results_go_dt
  }

# continuous versions ----------------------------------------------------------

#' Run GO enrichment with elimination with fgsea simple
#'
#' @description
#' This method takes the gene_ontology_data and a vector of gene level
#' statistics to perform fgsea (simple) leveraging ontological information.
#' It starts at the lowest levels of the ontology and tests if there is
#' significant enrichment for any GO terms. If the threshold of the p-value is
#' below the elimination threshold, the genes from this term will be removed
#' from all its ancestors. The function then proceeds to the next level of the
#' ontology and repeats the process.
#'
#' @param object The underlying class, see [bixverse::gene_ontology_data()].
#' @param stats Named numeric vector. The gene level statistic.
#' @param elim_threshold Float. Threshold from which p-value onwards the
#' elimination on the ancestors shall be conducted.
#' @param nperm Integer. Number of permutation tests. Defaults to `2000L`
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
#' @return data.table with enrichment results.
#'
#' @export
fgsea_simple_go_elim <- S7::new_generic(
  name = "fgsea_simple_go_elim",
  dispatch_args = "object",
  fun = function(
    object,
    stats,
    elim_threshold = 0.05,
    nperm = 2000L,
    gsea_params = params_gsea(max_size = 2000L),
    seed = 42L
  ) {
    S7::S7_dispatch()
  }
)

#' @export
#'
#' @import data.table
#' @importFrom magrittr `%>%`
#'
#' @method fgsea_simple_go_elim gene_ontology_data
S7::method(fgsea_simple_go_elim, gene_ontology_data) <-
  function(
    object,
    stats,
    elim_threshold = 0.05,
    nperm = 2000L,
    gsea_params = params_gsea(max_size = 2000L),
    seed = 42L
  ) {
    # Scope checks
    . <- leading_edge <- fdr <- NULL
    # Checks
    checkmate::assertClass(object, "bixverse::gene_ontology_data")
    checkmate::assertNumeric(stats, min.len = 3L, finite = TRUE)
    checkmate::assertNames(names(stats))
    checkmate::qassert(elim_threshold, "R+[0,1]")
    checkmate::qassert(nperm, "I1")
    checkmate::qassert(seed, "I1")
    assertGSEAParams(gsea_params)

    # Function body
    stats <- sort(stats, decreasing = TRUE)

    levels <- names(S7::prop(object, "levels"))
    go_info <- S7::prop(object, "go_info")

    rust_res <- rs_geom_elim_fgsea_simple(
      stats = stats,
      levels = levels,
      go_obj = object,
      gsea_params = gsea_params,
      elim_threshold = elim_threshold,
      iters = nperm,
      seed = seed
    )

    leading_edges <- mapply(
      "[",
      list(names(stats)),
      rust_res$leading_edge,
      SIMPLIFY = FALSE
    )

    res_dt <- data.table::setDT(rust_res[c(
      "go_id",
      "es",
      "nes",
      "size",
      "pvals"
    )]) %>%
      .[, leading_edge := leading_edges] %>%
      merge(., go_info, by = "go_id") %>%
      .[, fdr := p.adjust(pvals, method = "fdr")]

    return(res_dt)
  }


#' Run GO enrichment with elimination method over a continuous vectors
#'
#' @description
#' This method takes the gene_ontology_data and a vector of gene level
#' statistics to perform fgsea (multi-level) leveraging ontological information.
#' It starts at the lowest levels of the ontology and tests if there is
#' significant enrichment for any GO terms. If the threshold of the p-value is
#' below the elimination threshold, the genes from this term will be removed
#' from all its ancestors. The function then proceeds to the next level of the
#' ontology and repeats the process. Subsequently, it leverages the multi-level
#' method to estimate lower p-values for significant terms, see Korotkevich, et
#' al.
#'
#' @param object The underlying class, see [bixverse::gene_ontology_data()].
#' @param stats Named numeric vector. The gene level statistic.
#' @param elim_threshold Float. Threshold from which p-value onwards the
#' elimination on the ancestors shall be conducted.
#' @param nperm Integer. Number of permutation tests. Defaults to `2000L`
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
#' @return data.table with enrichment results.
#'
#' @references Korotkevich, et al., bioRxiv
#'
#' @export
fgsea_go_elim <- S7::new_generic(
  name = "fgsea_go_elim",
  dispatch_args = "object",
  fun = function(
    object,
    stats,
    elim_threshold = 0.05,
    nperm = 2000L,
    gsea_params = params_gsea(max_size = 2000L),
    seed = 42L
  ) {
    S7::S7_dispatch()
  }
)

#' @export
#'
#' @import data.table
#' @importFrom magrittr `%>%`
#'
#' @method fgsea_go_elim gene_ontology_data
S7::method(fgsea_go_elim, gene_ontology_data) <-
  function(
    object,
    stats,
    elim_threshold = 0.05,
    nperm = 2000L,
    gsea_params = params_gsea(max_size = 2000L),
    seed = 42L
  ) {
    # Scope checks
    . <- leading_edge <- fdr <- mode_fraction <- pvals <- log2err <- NULL
    # Checks
    checkmate::assertClass(object, "bixverse::gene_ontology_data")
    checkmate::assertNumeric(stats, min.len = 3L, finite = TRUE)
    checkmate::assertNames(names(stats))
    checkmate::qassert(elim_threshold, "R+[0,1]")
    checkmate::qassert(nperm, "I1")
    checkmate::qassert(seed, "I1")
    assertGSEAParams(gsea_params)

    # Function body
    stats <- sort(stats, decreasing = TRUE)

    levels <- names(S7::prop(object, "levels"))
    go_info <- S7::prop(object, "go_info")

    rust_res_simple <- rs_geom_elim_fgsea_simple(
      stats = stats,
      levels = levels,
      go_obj = object,
      gsea_params = gsea_params,
      elim_threshold = elim_threshold,
      iters = nperm,
      seed = seed
    ) %>%
      data.table::setDT()

    leading_edges <- mapply(
      "[",
      list(names(stats)),
      rust_res_simple$leading_edge,
      SIMPLIFY = FALSE
    )

    rust_res_simple[, `:=`(
      leading_edge = leading_edges,
      mode_fraction = data.table::fifelse(es >= 0, ge_zero, le_zero)
    )]

    rs_err_res <- with(
      gsea_params,
      rs_simple_and_multi_err(
        n_more_extreme = as.integer(rust_res_simple$n_more_extreme),
        nperm = nperm,
        sample_size = sample_size
      )
    )

    dt_simple_gsea <- rust_res_simple[
      rs_err_res$multi_err >= rs_err_res$simple_err
    ][,
      `:=`(
        log2err = 1 /
          log(2) *
          sqrt(trigamma(n_more_extreme + 1) - trigamma(nperm + 1))
      )
    ]

    dt_multi_level <- rust_res_simple[
      rs_err_res$multi_err < rs_err_res$simple_err
    ][, "denom_prob" := (mode_fraction + 1) / (nperm + 1)]

    rs_multi_level_res = with(
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

    dt_multi_level = dt_multi_level[, pvals := rs_multi_level_res$pvals] %>%
      data.table::as.data.table() %>%
      .[, pvals := pmin(1, pvals / denom_prob)] %>%
      .[,
        log2err := multilevel_error(
          pvals,
          sample_size = gsea_params$sample_size
        )
      ] %>%
      .[,
        log2err := data.table::fifelse(
          rs_multi_level_res$is_cp_ge_half,
          log2err,
          NA
        )
      ]

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
      merge(., go_info, by = "go_id") %>%
      data.table::setorder(pvals)

    return(all_results)
  }
