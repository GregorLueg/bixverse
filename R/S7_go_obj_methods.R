# Gene Set Enrichment versions ----

#' GSE on a single set of target genes with the GO elimination method.
#'
#' @description
#' This method takes the gene_ontology_data and a target gene set and performs an GSE
#' enrichment leveraging ontological information. It starts at the lowest levels of the ontology
#' and tests if there is significant enrichment for any GO terms. If the threshold of the p-value is below
#' the elimination threshold, the genes from this term will be removed from all its ancestors. The
#' function then proceeds to the next level of the ontology and repeats the process.
#'
#' @usage GSE_GO_elim_method(
#'  S7_obj,
#'  target_genes,
#'  minimum_overlap = 3L,
#'  fdr_threshold = 0.05,
#'  min_genes = NULL,
#'  .debug = FALSE
#' )
#'
#' @param S7_obj `gene_ontology_data` The underlying gene ontology data.
#' @param target_genes String. The target genes you wish to apply the GSEA over.
#' @param minimum_overlap Integer. Threshold for the minimal overlap.
#' @param fdr_threshold Float. Threshold for maximum FDR to include in the output.
#' @param elim_threshold Float. Threshold from which p-value onwards the elimination on the ancestors shall
#' be conducted.
#' @param min_genes Integer. Minimum number of genes that have to be included in the gene ontology term
#' @param .debug Boolean. Shall information from the Rust function be displayed. For debugging purposes.
#'
#' @return data.table with enrichment results.
#'
#' @export
GSE_GO_elim_method <- S7::new_generic("GSE_GO_elim_method", "S7_obj")

S7::method(GSE_GO_elim_method, gene_ontology_data) <-
  function(S7_obj,
           target_genes,
           minimum_overlap = 3L,
           fdr_threshold = 0.05,
           elim_threshold = 0.05,
           min_genes = NULL,
           .debug = FALSE) {
    # First check
    checkmate::assertClass(S7_obj, "BIXverse::gene_ontology_data")
    checkmate::qassert(target_genes, "S+")
    checkmate::qassert(fdr_threshold, "R+[0,1]")
    checkmate::qassert(elim_threshold, "R+[0,1]")
    checkmate::qassert(minimum_overlap, "I1")
    checkmate::qassert(min_genes, c("0", "I1"))
    checkmate::qassert(.debug, "B1")
    # Extract relevant data from the S7 object
    if (is.null(min_genes)) {
      min_genes <- S7::prop(S7_obj, "min_genes")
    }
    go_to_genes <- S7::prop(S7_obj, "go_to_genes")
    ancestry <- S7::prop(S7_obj, "ancestry")
    levels <- S7::prop(S7_obj, "levels")
    go_info <- S7::prop(S7_obj, "go_info")

    gene_universe_length <- length(unique(unlist(go_to_genes)))

    results_go <- rs_gse_geom_elim(
      target_genes = target_genes,
      go_to_genes = go_to_genes,
      ancestors = ancestry,
      levels = levels,
      gene_universe_length = gene_universe_length,
      min_genes = min_genes,
      elim_threshold = elim_threshold,
      debug = .debug
    )

    # Weird issue sometimes with pulling f64 back to Double and end up with negative p-values
    results_go_dt <- data.table(do.call(cbind, results_go[-1])) %>%
      .[, pvals := data.table::fifelse(pvals < 0, pvals * -1, pvals)] %>%
      .[, `:=`(go_id = results_go$go_ids,
               FDR = p.adjust(pvals, method = "BH"))] %>%
      data.table::setcolorder(.,
                              c(
                                "go_id",
                                "odds_ratios",
                                "pvals",
                                "FDR",
                                "hits",
                                "gene_set_lengths"
                              ))

    results_go_dt <- merge(go_info, results_go_dt, by = "go_id") %>%
      data.table::setorder(., pvals) %>%
      .[(FDR <= fdr_threshold) & (hits >= minimum_overlap)]

    results_go_dt
  }



#' GSE on a list of target genes with the GO elimination method.
#'
#' @description
#' This method takes the gene_ontology_data and a list of target gene sets and performs an GSE
#' enrichment leveraging ontological information. It starts at the lowest levels of the ontology
#' and tests if there is significant enrichment for any GO terms. If the threshold of the p-value is below
#' the elimination threshold, the genes from this term will be removed from all its ancestors. The
#' function then proceeds to the next level of the ontology and repeats the process. The class will
#' leverage Rust threading to parallelise the process.
#'
#' @usage GSE_GO_elim_method(
#'  S7_obj,
#'  target_gene_list,
#'  minimum_overlap = 3L,
#'  fdr_threshold = 0.05,
#'  min_genes = NULL,
#'  .debug = FALSE
#' )
#'
#' @param S7_obj `gene_ontology_data` The underlying gene ontology data.
#' @param target_gene_list List. The target genes list you wish to apply the GSEA over.
#' @param minimum_overlap Integer. Threshold for the minimal overlap.
#' @param fdr_threshold Float. Threshold for maximum FDR to include in the output.
#' @param elim_threshold Float. Threshold from which p-value onwards the elimination on the ancestors shall
#' be conducted.
#' @param min_genes Integer. Minimum number of genes that have to be included in the gene ontology term
#' @param .debug Boolean. Shall information from the Rust function be displayed. For debugging purposes.
#' Warning: should you run this command over a large list, you will have a large print output!
#'
#' @return data.table with enrichment results.
#'
#' @export
GSE_GO_elim_method_list <- S7::new_generic("GSE_GO_elim_method_list", "S7_obj")

S7::method(GSE_GO_elim_method_list, gene_ontology_data) <-
  function(S7_obj,
           target_gene_list,
           minimum_overlap = 3L,
           fdr_threshold = 0.05,
           elim_threshold = 0.05,
           min_genes = NULL,
           .debug = FALSE) {
    # First check
    checkmate::assertClass(S7_obj, "BIXverse::gene_ontology_data")
    checkmate::assertList(target_gene_list, types = "character")
    checkmate::qassert(fdr_threshold, "R+[0,1]")
    checkmate::qassert(elim_threshold, "R+[0,1]")
    checkmate::qassert(minimum_overlap, "I1")
    checkmate::qassert(min_genes, c("0", "I1"))
    checkmate::qassert(.debug, "B1")
    # Extract relevant data from the S7 object
    if (is.null(min_genes)) {
      min_genes <- S7::prop(S7_obj, "min_genes")
    }
    go_to_genes <- S7::prop(S7_obj, "go_to_genes")
    ancestry <- S7::prop(S7_obj, "ancestry")
    levels <- S7::prop(S7_obj, "levels")
    go_info <- S7::prop(S7_obj, "go_info")

    gene_universe_length <- length(unique(unlist(go_to_genes)))

    results_go <- rs_gse_geom_elim_list(
      target_genes = target_gene_list,
      go_to_genes = go_to_genes,
      ancestors = ancestry,
      levels = levels,
      gene_universe_length = gene_universe_length,
      min_genes = min_genes,
      elim_threshold = elim_threshold,
      debug = .debug
    )

    target_set_names <- unlist(purrr::map2(names(target_gene_list), results_go$no_test, ~ {
      rep(.x, each = .y)
    }))

    results_go_dt <- data.table(do.call(cbind, results_go[c("pvals", "odds_ratios", "hits", "gene_set_lengths")])) %>%
      .[, pvals := data.table::fifelse(pvals < 0, pvals * -1, pvals)] %>%
      .[, `:=`(go_id = results_go$go_ids, target_set_name = target_set_names)] %>%
      .[, FDR := p.adjust(pvals, method = "BH"), by = target_set_name]

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
          "FDR",
          "hits",
          "gene_set_lengths"
        )
      ) %>%
      .[(FDR <= fdr_threshold) & (hits >= minimum_overlap)]

    results_go_dt
  }
