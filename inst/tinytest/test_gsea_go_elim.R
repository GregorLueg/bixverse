# (fast) gsea test for gene ontology elimination -------------------------------

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
  gsea_param = 1.0,
  elim_threshold = 0.00001, # Threshold is so low, it cannot be passed
  min_size = 3,
  max_size = 250,
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

### wit elimination ------------------------------------------------------------

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
  gsea_param = 1.0,
  elim_threshold = 0.95, # This WILL be passed
  min_size = 3,
  max_size = 250,
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
