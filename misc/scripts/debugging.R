object = object
target_gene_list = go_target_genes
minimum_overlap = 0L
fdr_threshold = 1
elim_threshold = 0.95


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

cols_to_select <- c("pvals", "fdr", "odds_ratios", "hits", "gene_set_lengths")

results_go_dt <- data.table(do.call(cbind, rs_results_go[cols_to_select])) %>%
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
