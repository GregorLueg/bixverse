# gsea test for gene ontology elimination (simple) -----------------------------

library(magrittr)

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

toy_go_data <- data.table::data.table(
  go_id = sprintf("go_%i", 1:6),
  go_name = sprintf("go_name_%s", letters[1:6]),
  ancestors = list(
    c("go_1"),
    c("go_1", "go_2"),
    c("go_1", "go_2", "go_3"),
    c("go_4"),
    c("go_4", "go_5"),
    c("go_4", "go_5", "go_6")
  ),
  ensembl_id = list(
    pathway_pos_l3,
    pathway_pos_l2,
    pathway_pos_l1,
    pathway_neg_l3,
    pathway_neg_l2,
    pathway_neg_l1
  ),
  depth = c(1, 2, 3, 1, 2, 3)
) %>%
  data.table::setorder(-depth)

## test rust implementations ---------------------------------------------------

object <- gene_ontology_data(toy_go_data, min_genes = 3L)

levels <- names(S7::prop(object, "levels"))

### without elimination --------------------------------------------------------

expected_sizes_wo_elim <- c(35, 40, 64, 97, 100, 125)
expected_nes_wo_elim <- c(
  2.970346,
  -3.515364,
  3.472927,
  -3.759593,
  4.060339,
  -4.188695
)

rust_res_wo_elim <- rs_geom_elim_fgsea_simple(
  stats = stats,
  levels = levels,
  go_obj = object,
  gsea_params = params_gsea(min_size = 3L, max_size = 250L),
  elim_threshold = 0.00001, # Threshold is so low, it cannot be passed
  iters = 100,
  seed = 10101,
  debug = FALSE
)

expect_equal(
  current = rust_res_wo_elim$size,
  target = expected_sizes_wo_elim,
  info = paste(
    "fgsea with go elim (rust function): sizes in no elimination case"
  ),
  tolerance = 10e-6
)

expect_equal(
  current = rust_res_wo_elim$nes,
  target = expected_nes_wo_elim,
  info = paste(
    "fgsea with go elim (rust function): nes in no elimination case"
  ),
  tolerance = 10e-6
)

### with elimination -----------------------------------------------------------

# The elimination should yield smaller tested gs sizes and worse NES
# for later gene sets
expected_sizes_with_elim <- c(35, 40, 39, 62, 50, 46)
expected_nes_with_elim <- c(
  2.970346,
  -3.515364,
  2.689328,
  -3.014979,
  3.248230,
  -3.340476
)

rust_res_with_elim <- rs_geom_elim_fgsea_simple(
  stats = stats,
  levels = levels,
  go_obj = object,
  gsea_params = params_gsea(min_size = 3L, max_size = 250L),
  elim_threshold = 0.95, # This WILL be passed
  iters = 100,
  seed = 10101,
  debug = FALSE
)

expect_equal(
  current = rust_res_with_elim$size,
  target = expected_sizes_with_elim,
  info = paste(
    "fgsea with go elim (rust function): sizes in with elimination case"
  ),
  tolerance = 10e-6
)

expect_equal(
  current = rust_res_with_elim$nes,
  target = expected_nes_with_elim,
  info = paste(
    "fgsea with go elim (rust function): nes in with elimination case"
  ),
  tolerance = 10e-6
)

## test R into Rust implementation ---------------------------------------------

### without elimination --------------------------------------------------------

# order is different due to merging
r_expected_size_without_elim <- c(100, 64, 35, 125, 97, 40)
r_expected_nes_without_elim <- c(
  4.060339,
  3.472927,
  2.970346,
  -4.188695,
  -3.759593,
  -3.515364
)

r_results_without_elim <- fgsea_simple_go_elim(
  object = object,
  stats = stats,
  nperm = 100L,
  seed = 10101L,
  elim_threshold = 0.0001
)

expect_equal(
  current = r_results_without_elim$size,
  target = r_expected_size_without_elim,
  info = paste(
    "fgsea with go elim (R function): sizes in without elimination case"
  ),
  tolerance = 10e-6
)

expect_equal(
  current = r_results_without_elim$nes,
  target = r_expected_nes_without_elim,
  info = paste(
    "fgsea with go elim (R function): nes in with elimination case"
  ),
  tolerance = 10e-6
)

### with elimination -----------------------------------------------------------

r_expected_size_with_elim <- c(50, 39, 35, 46, 62, 40)
r_expected_nes_with_elim <- c(
  3.248230,
  2.689328,
  2.970346,
  -3.340476,
  -3.014979,
  -3.515364
)

r_results_with_elim <- fgsea_simple_go_elim(
  object = object,
  stats = stats,
  nperm = 100L,
  seed = 10101L,
  elim_threshold = 0.05
)

expect_equal(
  current = r_results_with_elim$size,
  target = r_expected_size_with_elim,
  info = paste(
    "fgsea with go elim (R function): sizes in with elimination case"
  ),
  tolerance = 10e-6
)

expect_equal(
  current = r_results_with_elim$nes,
  target = r_expected_nes_with_elim,
  info = paste(
    "fgsea with go elim (R function): nes in with elimination case"
  ),
  tolerance = 10e-6
)

# gsea test for gene ontology elimination (multi level) ------------------------

## data ------------------------------------------------------------------------

# data set with non-significant genes to test the
# multi-level p-value calculations

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

## test R implementation -------------------------------------------------------

### with elimination -----------------------------------------------------------

# these are super significant due to the sampling...

expected_sizes <- c(35, 48, 61, 44, 40, 38, 40, 50)
expected_err <- c(NA, NA, NA, NA, NA, 1.00087305, 0.06299909, 0.03458907)

r_results <- fgsea_go_elim(
  object = object,
  stats = stats,
  elim_threshold = 0.95
)

expect_equal(
  current = r_results$size,
  target = expected_sizes,
  info = paste(
    "fgsea (multi level) with go elim (R function):",
    "sizes in with elimination case"
  ),
  tolerance = 10e-6
)

expect_equal(
  current = r_results$log2err,
  target = expected_err,
  info = paste(
    "fgsea (multi level) with go elim (R function):",
    "log2err in with elimination case"
  ),
  tolerance = 10e-6
)

### without elimination --------------------------------------------------------

expected_sizes <- c(35, 100, 64, 97, 40, 125, 40, 50)
expected_err <- c(NA, NA, NA, NA, NA, NA, 0.06299909, 0.03458907)

r_results <- fgsea_go_elim(
  object = object,
  stats = stats,
  elim_threshold = 0.00001
)

expect_equal(
  current = r_results$size,
  target = expected_sizes,
  info = paste(
    "fgsea (multi level) with go elim (R function):",
    "sizes in without elimination case"
  ),
  tolerance = 10e-6
)

expect_equal(
  current = r_results$log2err,
  target = expected_err,
  info = paste(
    "fgsea (multi level) with go elim (R function):",
    "log2err in without elimination case"
  ),
  tolerance = 10e-6
)
