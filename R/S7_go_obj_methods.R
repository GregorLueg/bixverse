#' GSE on a single set of target genes with the GO elimination method
#'
#' @description
#' ...
#'
#' @param S7_obj `gene_ontology_data` class.
#' @param target_genes String. ...
#' @param minimum_overlap Integer. ...
#' @param fdr_threshold Float. ...
#' @param elim_threshold Float. ...
#' @param min_genes Integer. ...
#' @param .debug Boolean. ...
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
    # checkmate::assertClass(S7_obj, "gene_ontology_data")
    checkmate::qassert(target_genes, "S+")
    checkmate::qassert(fdr_threshold, "R+[0,1]")
    checkmate::qassert(elim_threshold, "R+[0,1]")
    checkmate::qassert(minimum_overlap, "I1")
    checkmate::qassert(min_genes, c("0", "I1"))
    checkmate::qassert(.debug, "B1")
    # Extract relevant data from the S7 object
    if(is.null(min_genes)) {
      min_genes = S7::prop(S7_obj, "min_genes")
    }
    go_to_genes = S7::prop(S7_obj, "go_to_genes")
    ancestry = S7::prop(S7_obj, "ancestry")
    levels = S7::prop(S7_obj, "levels")
    go_info = S7::prop(S7_obj, "go_info")

    gene_universe_length = length(unique(unlist(go_to_genes)))

    results_go = rs_gse_geom_elim(
      target_genes = target_genes,
      go_to_genes = go_to_genes,
      ancestors = ancestry,
      levels = levels,
      gene_universe_length = gene_universe_length,
      min_genes = min_genes,
      elim_threshold = elim_threshold,
      debug = .debug
    )

    results_go_dt = data.table(do.call(cbind, results_go[-1])) %>%
      .[, `:=`(go_id = results_go$go_ids,
               FDR = p.adjust(pvals, method = 'BH'))] %>%
      data.table::setcolorder(
        .,
        c(
          'go_id',
          'odds_ratios',
          'pvals',
          'FDR',
          'hits',
          'gene_set_lengths'
        )
      )

    results_go_dt = merge(go_info,
                          results_go_dt,
                          by = 'go_id') %>%
      data.table::setorder(., pvals) %>%
      .[(FDR <= fdr_threshold) & (hits >= minimum_overlap)]

    results_go_dt
  }
