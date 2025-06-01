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


# Original fgsea multi level code -----

set.seed(123)

stat_size <- 20000

stats <- setNames(
  sort(rnorm(stat_size), decreasing = TRUE),
  paste0("gene", 1:stat_size)
)

number_gene_sets <- 500
min_size <- 50
max_size <- 250

pathway_random <- purrr::map(
  seq_len(number_gene_sets),
  ~ {
    sample(names(stats), sample(min_size:max_size, 1))
  }
)

names(pathway_random) <- paste0("pathway", 1:number_gene_sets)


pathways = pathway_random
stats = stats
sampleSize = 121
minSize = 1
maxSize = length(stats) - 1
eps = 1e-50
scoreType = "std"
nproc = 0
gseaParam = 1
BPPARAM = NULL
nPermSimple = 1000
absEps = NULL

pp <- fgsea:::preparePathwaysAndStats(
  pathways,
  stats,
  minSize,
  maxSize,
  gseaParam,
  scoreType
)

pathwaysFiltered <- pp$filtered
pathwaysSizes <- pp$sizes
stats <- pp$stats
m <- length(pathwaysFiltered)

minSize <- max(minSize, 1)
eps <- max(0, min(1, eps))

if (sampleSize %% 2 == 0) {
  sampleSize <- sampleSize + 1
}

gseaStatRes <- do.call(
  rbind,
  lapply(
    pathwaysFiltered,
    fgsea::calcGseaStat,
    stats = stats,
    returnLeadingEdge = TRUE,
    scoreType = scoreType
  )
)

leadingEdges <- mapply(
  "[",
  list(names(stats)),
  gseaStatRes[, "leadingEdge"],
  SIMPLIFY = FALSE
)
pathwayScores <- unlist(gseaStatRes[, "res"])

seeds <- sample.int(10^9, 1)

simpleFgseaRes <- fgsea:::fgseaSimpleImpl(
  pathwayScores = pathwayScores,
  pathwaysSizes = pathwaysSizes,
  pathwaysFiltered = pathwaysFiltered,
  leadingEdges = leadingEdges,
  permPerProc = nPermSimple,
  seeds = seeds,
  toKeepLength = m,
  stats = stats,
  BPPARAM = BiocParallel::SerialParam(),
  scoreType = scoreType
)

switch(
  scoreType,
  std = simpleFgseaRes[, modeFraction := ifelse(ES >= 0, nGeZero, nLeZero)],
  pos = simpleFgseaRes[, modeFraction := nGeZero],
  neg = simpleFgseaRes[, modeFraction := nLeZero]
)

simpleFgseaRes[, leZeroMean := NULL]
simpleFgseaRes[, geZeroMean := NULL]
simpleFgseaRes[, nLeEs := NULL]
simpleFgseaRes[, nGeEs := NULL]
simpleFgseaRes[, nLeZero := NULL]
simpleFgseaRes[, nGeZero := NULL]

rextendr::document()

rs_err_res <- rs_simple_and_multi_err(
  n_more_extreme = as.integer(simpleFgseaRes$nMoreExtreme),
  nperm = 1000L,
  sample_size = 121L
)


leftBorder <- log2(qbeta(
  0.025,
  shape1 = simpleFgseaRes$nMoreExtreme,
  shape2 = nPermSimple - simpleFgseaRes$nMoreExtreme + 1
))

rightBorder <- log2(qbeta(
  1 - 0.025,
  shape1 = simpleFgseaRes$nMoreExtreme + 1,
  shape2 = nPermSimple - simpleFgseaRes$nMoreExtreme
))

crudeEstimator <- log2((simpleFgseaRes$nMoreExtreme + 1) / (nPermSimple + 1))

simpleError <- 0.5 *
  pmax(crudeEstimator - leftBorder, rightBorder - crudeEstimator)


multError <- sapply(
  (simpleFgseaRes$nMoreExtreme + 1) / (nPermSimple + 1),
  fgsea::multilevelError,
  sampleSize
)

plot(
  rs_err_res$simple_err,
  simpleError
)

plot(
  rs_err_res$multi_err,
  multError
)

simpleError

# We do not bother taking these forward as the pvals >= 0.05
dtSimpleFgsea <- simpleFgseaRes[multError >= simpleError]
dtSimpleFgsea[,
  log2err := 1 /
    log(2) *
    sqrt(trigamma(nMoreExtreme + 1) - trigamma(nPermSimple + 1))
]
dtSimpleFgsea[, modeFraction := NULL]


dtMultilevel <- simpleFgseaRes[multError < simpleError]
dtMultilevel[, "denomProb" := (modeFraction + 1) / (nPermSimple + 1)]

multilevelPathwaysList <- split(dtMultilevel, by = "size")

indxs <- sample(1:length(multilevelPathwaysList))
multilevelPathwaysList <- multilevelPathwaysList[indxs]

seed = sample.int(1e9, size = 1)

sign <- if (scoreType %in% c("pos", "neg")) TRUE else FALSE

# This is where the magic is happening...

pathways = unlist(purrr::map(
  multilevelPathwaysList,
  ~ {
    .x$pathway
  }
))

multilevelPathwaysList

cpp.res <- fgsea:::multilevelImpl(
  multilevelPathwaysList,
  stats,
  121,
  seed,
  eps,
  sign = sign,
  BPPARAM = BiocParallel::SerialParam()
)

library(data.table)
library(magrittr)

cpp_res <- purrr::imap_dfr(cpp.res, \(val, name) {
  val <- data.table::as.data.table(val)
  val[, size := name]
}) %>%
  as.data.table()

length(pathways)

cpp_res[, pathway := pathways]

rextendr::document()


rs_res <- purrr::map(multilevelPathwaysList, \(x) {
  as.data.table(rs_calc_multi_level(
    es = x[, ES],
    stats = stats,
    pathway_size = unique(x[, size]),
    sample_size = 501,
    seed = 10101L,
    eps = eps,
    sign = sign
  ))
})

rs_res_dt = rbindlist(rs_res)[, pathway := pathways]

rs_res_dt

cpp_res

plot(cpp_res$cppMPval, rs_res_dt$pvals, ylim = c(0, 0.05), xlim = c(0, 0.05))

cor(cpp_res$cppMPval, rs_res_dt$pvals, method = 'spearman')

names(rs_res)

rs_res[["106"]]

cpp.res[["106"]]


rextendr::document()

tmp <- qs2::qs_read("~/Desktop/test.qs2")

purrr::iwalk(tmp, \(val, name) {
  assign(name, val, envir = .GlobalEnv)
})

rs_res = rs_calc_multi_level(
  es = x[, ES],
  stats = stats,
  pathway_size = unique(x[, size]),
  sample_size = sampleSize,
  seed = seed,
  eps = eps,
  sign = sign
)

## Implement a full function ---------------------------------------------------

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

library(zeallot)

stats = stats
pathways = pathway_random
nperm = 2000L
gsea_params = params_gsea()
seed = 123L
sign_ = TRUE

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

# Calculate the multi-level version
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

rextendr::document()

rs_res = with(
  gsea_params,
  rs_calc_multi_level(
    stats = stats,
    es = dt_multi_level$es,
    pathway_size = as.integer(dt_multi_level$size),
    sample_size = sample_size,
    seed = seed,
    eps = eps,
    sign = sign_
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
    denom_prob = NULL
  )]
