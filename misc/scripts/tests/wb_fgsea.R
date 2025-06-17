# Rust implementation of fgsea ----

set.seed(123)

stat_size <- 20000

stats <- setNames(
  sort(rnorm(stat_size), decreasing = TRUE),
  paste0("gene", 1:stat_size)
)

number_gene_sets <- 5000
min_size <- 50
max_size <- 250

pathway_random <- purrr::map(
  seq_len(number_gene_sets),
  ~ {
    sample(names(stats), sample(min_size:max_size, 1))
  }
)

names(pathway_random) <- paste0("pathway", 1:number_gene_sets)

devtools::load_all()

tictoc::tic()
results_traditional <- calc_gsea_traditional(
  stats = stats,
  pathways = pathway_random
)
tictoc::toc()

tictoc::tic()
results_simple_fgsea <- calc_fgsea_simple(
  stats = stats,
  pathways = pathway_random
)
tictoc::toc()


tictoc::tic()
fgsea_scores_original <- fgsea::fgseaSimple(
  pathways = pathway_random,
  stats = stats,
  nperm = 2000L
)
tictoc::toc()

tictoc::tic()
results_ml_fgsea <- calc_fgsea(
  stats = stats,
  pathways = pathway_random
)
tictoc::toc()


tictoc::tic()
fgsea_scores_ml <- fgsea::fgsea(
  stats = stats,
  pathways = pathway_random
)
tictoc::toc()

# multi-level fgsea go elimination ---------------------------------------------

## data ------------------------------------------------------------------------

set.seed(123)

stat_size <- 1000

stats <- setNames(
  sort(rnorm(stat_size), decreasing = TRUE),
  sprintf("gene_%03d", 1:stat_size)
)

pathway_pos_l3 <- sort(sample(names(stats)[1:150], 100))
pathway_pos_l2 <- unique(c(pathway_pos_l3[1:50], sample(names(stats), 15)))
pathway_pos_l1 <- unique(c(pathway_pos_l2[1:25], sample(names(stats), 10)))

pathway_neg_l3 <- sort(
  sample(names(stats)[851:1000], 125),
  decreasing = TRUE
)
pathway_neg_l2 <- unique(c(pathway_neg_l3[1:75], sample(names(stats), 25)))
pathway_neg_l1 <- unique(c(pathway_neg_l2[1:35], sample(names(stats), 5)))

pathway_random_l1.1 <- sample(names(stats), 40)
pathway_random_l1.2 <- sample(names(stats), 50)

toy_go_data <- data.table::data.table(
  go_id = sprintf("go_%i", 1:8),
  go_name = sprintf("go_name_%s", letters[1:8]),
  ancestors = list(
    c("go_1"),
    c("go_1", "go_2"),
    c("go_1", "go_2", "go_3"),
    c("go_1", "go_2", "go_4"),
    c("go_5"),
    c("go_5", "go_6"),
    c("go_5", "go_6", "go_7"),
    c("go_5", "go_6", "go_8")
  ),
  ensembl_id = list(
    pathway_pos_l3,
    pathway_pos_l2,
    pathway_pos_l1,
    pathway_random_l1.1,
    pathway_neg_l3,
    pathway_neg_l2,
    pathway_neg_l1,
    pathway_random_l1.2
  ),
  depth = c(1, 2, 3, 3, 1, 2, 3, 3)
) %>%
  data.table::setorder(-depth)

object <- gene_ontology_data(toy_go_data, min_genes = 3L)

levels <- names(S7::prop(object, "levels"))

## workbench -------------------------------------------------------------------

devtools::load_all()

results <- fgsea_go_elim(
  object = object,
  stats = stats,
  elim_threshold = 0.95
)

gsea_params <- params_gsea(min_size = 3L, max_size = 250L)
nperm = 1000
seed = 10101L

rust_res_simple <- rs_geom_elim_fgsea_simple(
  stats = stats,
  levels = levels,
  go_obj = object,
  gsea_params = gsea_params,
  elim_threshold = 0.95,
  iters = nperm,
  seed = seed,
  debug = FALSE
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
    log2err := multilevel_error(pvals, sample_size = gsea_params$sample_size)
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
  data.table::setorder(pvals) %>%
  merge(., go_info, by = "go_id")
