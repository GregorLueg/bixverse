# library ----------------------------------------------------------------------

library(magrittr)
library(zeallot)

# test ontology semantic similarity --------------------------------------------

## expected values -------------------------------------------------------------

test_onto <- data.table::data.table(
  parent = c("a", "b", "b", "b", "c"),
  child = c("b", "c", "d", "e", "f")
)

expected_ancestors <- list(
  a = "a",
  b = c("b", "a"),
  c = c("c", "b", "a"),
  d = c("d", "b", "a"),
  e = c("e", "b", "a"),
  f = c("f", "c", "b", "a")
)

expected_descendants <- list(
  "a" = c("b", "c", "d", "e", "f", "a"),
  "b" = c("b", "c", "d", "e", "f"),
  "c" = c("c", "f"),
  "d" = c("d"),
  "e" = c("e"),
  "f" = c("f")
)

expected_ic_data <- list(
  a = 0,
  b = 0.1823216,
  c = 1.098612,
  d = 1.791759,
  e = 1.791759,
  f = 1.791759
)

expected_resnik <- c(
  0,
  0,
  0,
  0,
  0,
  0.1823216,
  0.1823216,
  0.1823216,
  0.1823216,
  0.1823216,
  0.1823216,
  1.0986123,
  0.1823216,
  0.1823216,
  0.1823216
)

expected_lin <- c(
  0,
  0,
  0,
  0,
  0,
  0.2846697,
  0.1847154,
  0.1847154,
  0.1847154,
  0.1261579,
  0.1261579,
  0.7601875,
  0.1017556,
  0.1017556,
  0.1017556
)

expected_combined <- c(
  0,
  0,
  0,
  0,
  0,
  0.1932127,
  0.1432355,
  0.1432355,
  0.1432355,
  0.1139567,
  0.1139567,
  0.6866674,
  0.1017556,
  0.1017556,
  0.1017556
)

expected_data_table <- data.table::data.table(
  feature_a = c(
    "a",
    "a",
    "a",
    "a",
    "a",
    "b",
    "b",
    "b",
    "b",
    "c",
    "c",
    "c",
    "d",
    "d",
    "e"
  ),
  feature_b = c(
    "b",
    "c",
    "d",
    "e",
    "f",
    "c",
    "d",
    "e",
    "f",
    "d",
    "e",
    "f",
    "e",
    "f",
    "f"
  ),
  sim = expected_resnik
)

## separate functions ----------------------------------------------------------

### ancestry ------------------------------------------------------------------

c(ancestors, descendants) %<-% get_ontology_ancestry(test_onto)

expect_equivalent(
  current = ancestors[names(expected_ancestors)],
  target = expected_ancestors,
  info = "Ontology similarity ancestors test"
)

expect_equivalent(
  current = descendants[names(expected_descendants)],
  target = expected_descendants,
  info = "Ontology similarity descendants test"
)

### information content --------------------------------------------------------

ic_data <- calculate_information_content(descendants)

expect_equivalent(
  current = ic_data[names(expected_ic_data)],
  target = expected_ic_data,
  info = "Ontology similarity test for information content.",
  tolerance = 1e-6
)

### semantic similarity --------------------------------------------------------

resnik <- calculate_semantic_sim(
  similarity_type = "resnik",
  terms = sort(names(ancestors)),
  ancestor_list = ancestors,
  ic_list = ic_data
)

lin <- calculate_semantic_sim(
  similarity_type = "lin",
  terms = sort(names(ancestors)),
  ancestor_list = ancestors,
  ic_list = ic_data
)

combined <- calculate_semantic_sim(
  similarity_type = "combined",
  terms = sort(names(ancestors)),
  ancestor_list = ancestors,
  ic_list = ic_data
)

expect_equivalent(
  current = rs_dense_to_upper_triangle(resnik, 1),
  target = expected_resnik,
  info = "Ontology similarity test for semantic semilarity (Resnik).",
  tolerance = 1e-6
)

expect_equivalent(
  current = rs_dense_to_upper_triangle(lin, 1),
  target = expected_lin,
  info = "Ontology similarity test for semantic semilarity (Lin).",
  tolerance = 1e-6
)

expect_equivalent(
  current = rs_dense_to_upper_triangle(combined, 1),
  target = expected_combined,
  info = "Ontology similarity test for semantic semilarity (combined type).",
  tolerance = 1e-6
)

## class -----------------------------------------------------------------------

test_class <- ontology(test_onto, .verbose = FALSE)

expect_warning(
  current = calculate_semantic_sim_onto(test_class, sim_type = "resnik"),
  info = paste("Ontology similarity - warning working")
)

test_class <- pre_process_sim_onto(test_class, .verbose = FALSE)

test_class <- calculate_semantic_sim_onto(object = test_class, .verbose = FALSE)

matrix_result <- get_sim_matrix(test_class, .verbose = FALSE)

dt_result <- get_sim_matrix(test_class, as_data_table = TRUE, .verbose = FALSE)


expect_equivalent(
  current = rs_dense_to_upper_triangle(matrix_result, 1),
  target = expected_resnik,
  info = paste(
    "Ontology class semantic similarity - matrix version"
  ),
  tolerance = 1e-6
)

expect_equal(
  current = dt_result,
  target = expected_data_table,
  info = paste(
    "Ontology class semantic similarity - data.table version"
  ),
  tolerance = 1e-6
)

# test ontology wang similarity ------------------------------------------------

## expected values -------------------------------------------------------------

test_onto <- data.table::data.table(
  parent = c("a", "b", "b", "b", "c", "g"),
  child = c("b", "c", "d", "e", "f", "h")
)

expected_sims_w08 <- c(
  0.4807122,
  0,
  0.6212121,
  0.3825911,
  0.8145401,
  0,
  0.4807122,
  0,
  0.7641509,
  0.4767442,
  0.5901639,
  0,
  0.5901639,
  0,
  0,
  0,
  0.6428571,
  0,
  0.6428571,
  0.7641509,
  0,
  0.7641509,
  0.4767442,
  0,
  0.4767442,
  0,
  0.5901639,
  0
)

expected_sims_w06 <- c(
  0.3713733,
  0,
  0.5762712,
  0.3828715,
  0.7582205,
  0,
  0.3713733,
  0,
  0.7191011,
  0.4594595,
  0.4897959,
  0,
  0.4897959,
  0,
  0,
  0,
  0.6153846,
  0,
  0.6153846,
  0.7191011,
  0,
  0.7191011,
  0.4594595,
  0,
  0.4594595,
  0,
  0.4897959,
  0
)

expected_dt_w08 <- data.table::data.table(
  feature_a = c(
    "f",
    "f",
    "f",
    "f",
    "f",
    "f",
    "f",
    "d",
    "d",
    "d",
    "d",
    "d",
    "d",
    "g",
    "g",
    "g",
    "g",
    "g",
    "b",
    "b",
    "b",
    "b",
    "a",
    "a",
    "a",
    "c",
    "c",
    "h"
  ),
  feature_b = c(
    "d",
    "g",
    "b",
    "a",
    "c",
    "h",
    "e",
    "g",
    "b",
    "a",
    "c",
    "h",
    "e",
    "b",
    "a",
    "c",
    "h",
    "e",
    "a",
    "c",
    "h",
    "e",
    "c",
    "h",
    "e",
    "h",
    "e",
    "e"
  ),
  sim = expected_sims_w08
)

expected_critval_w08 <- 0.7641509

expected_critval_w06 <- 0.7191011

## functions -------------------------------------------------------------------

results_w08 <- calculate_wang_sim(test_onto, w = 0.8)

results_w06 <- calculate_wang_sim(test_onto, w = 0.6)

critval_w08 <- calculate_critical_value(results_w08, alpha = 0.1)

critval_w06 <- calculate_critical_value(results_w06, alpha = 0.1)

expect_equivalent(
  current = rs_dense_to_upper_triangle(results_w08, 1L),
  target = expected_sims_w08,
  info = "Wang similarity w = 0.8 test",
  tolerance = 1e-6
)

expect_equivalent(
  current = rs_dense_to_upper_triangle(results_w06, 1L),
  target = expected_sims_w06,
  info = "Wang similarity w = 0.6 test",
  tolerance = 1e-6
)

expect_equivalent(
  current = critval_w08,
  target = expected_critval_w08,
  info = "Wang similarity w = 0.8 critical value",
  tolerance = 1e-6
)

expect_equivalent(
  current = critval_w06,
  target = expected_critval_w06,
  info = "Wang similarity w = 0.6 critical value",
  tolerance = 1e-6
)

## class -----------------------------------------------------------------------

test_class <- ontology(test_onto, .verbose = FALSE)

test_class <- calculate_wang_sim_onto(test_class)

matrix_res <- get_sim_matrix(
  test_class,
  as_data_table = FALSE,
  .verbose = FALSE
)

dt_res <- get_sim_matrix(
  test_class,
  as_data_table = TRUE,
  .verbose = FALSE
)

critval_class <- calculate_critical_value(test_class, alpha = 0.1)

expect_equivalent(
  current = rs_dense_to_upper_triangle(matrix_res, 1L),
  target = expected_sims_w08,
  info = "Ontology class wang similarity - matrix version",
  tolerance = 1e-6
)

expect_equal(
  current = dt_res,
  target = expected_dt_w08,
  info = paste(
    "Ontology class wang similarity - data.table version"
  ),
  tolerance = 1e-6
)

expect_equivalent(
  current = critval_class,
  target = expected_critval_w08,
  info = "Ontology class wang similarity - critical value",
  tolerance = 1e-6
)
