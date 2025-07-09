library("fgsea")

data(examplePathways)
data(exampleRanks)

stats = exampleRanks
pathways = examplePathways
gsea_params = params_gsea(min_size = 15L)
nperm = 1000L
seed = 123L

tictoc::tic()
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
tictoc::toc()

# 0.029 sec elapsed

tictoc::tic()
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
tictoc::toc()

# 0.008 sec elapsed

tictoc::tic()
leading_edges <- mapply(
  "[",
  list(names(stats)),
  gsea_stat_res[, "leading_edge"],
  SIMPLIFY = FALSE
)
tictoc::toc()

# 0.002 sec elapsed

pathway_scores <- unlist(gsea_stat_res[, "es"])

tictoc::tic()
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
tictoc::toc()

# 0.046 sec elapsed

# Calculations for the multi-level version
tictoc::tic()
rs_err_res <- with(
  gsea_params,
  rs_simple_and_multi_err(
    n_more_extreme = as.integer(permutations_res_simple$n_more_extreme),
    nperm = nperm,
    sample_size = sample_size
  )
)
tictoc::toc()

# 0.002 sec elapsed

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

tictoc::tic()
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
tictoc::toc()

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
