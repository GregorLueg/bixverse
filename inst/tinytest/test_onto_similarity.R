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
  0.1017556,
  0.1017556,
  0.1017556,
  0.1017556,
  0.1017556,
  0.1017556,
  0.6131472,
  0.1017556,
  0.1017556,
  0.1017556
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

expected_subset <- data.table::data.table(
  term1 = c("a", "a", "c"),
  term2 = c("c", "f", "f"),
  sims = c(0, 0, 0.6131472)
)

expected_filtered <- data.table::data.table(
  t1 = "c",
  t2 = "f",
  sim = 0.6131472
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

#### full matrices -------------------------------------------------------------

resnik <- calculate_semantic_sim_mat(
  similarity_type = "resnik",
  ancestor_list = ancestors,
  ic_list = ic_data
)

lin <- calculate_semantic_sim_mat(
  similarity_type = "lin",
  ancestor_list = ancestors,
  ic_list = ic_data
)

combined <- calculate_semantic_sim_mat(
  similarity_type = "combined",
  ancestor_list = ancestors,
  ic_list = ic_data
)

expect_equal(
  current = rs_dense_to_upper_triangle(resnik, 1),
  target = expected_resnik,
  info = "Ontology similarity test for semantic semilarity (Resnik).",
  tolerance = 1e-6
)

expect_equal(
  current = rs_dense_to_upper_triangle(lin, 1),
  target = expected_lin,
  info = "Ontology similarity test for semantic semilarity (Lin).",
  tolerance = 1e-6
)

expect_equal(
  current = rs_dense_to_upper_triangle(combined, 1),
  target = expected_combined,
  info = "Ontology similarity test for semantic semilarity (combined type).",
  tolerance = 1e-6
)

#### sub sets ------------------------------------------------------------------

resnik_subset <- calculate_semantic_sim(
  terms = c("a", "c", "f"),
  similarity_type = "resnik",
  ancestor_list = ancestors,
  ic_list = ic_data
)

expect_equal(
  current = resnik_subset,
  target = expected_subset,
  info = "Ontology similarity test for semantic semilarity (Resnik - subset).",
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

expect_equal(
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

### filtering ------------------------------------------------------------------

test_class <- filter_similarities(
  object = test_class,
  alpha = 0.01,
  .verbose = FALSE
)

filtered_results <- get_results(test_class)

expect_equal(
  current = filtered_results,
  target = expected_filtered,
  info = paste(
    "Ontology class semantic similarity - data.table version"
  ),
  tolerance = 1e-6
)

# test ontology wang similarity ------------------------------------------------

## expected values -------------------------------------------------------------

test_onto_wang <- data.table::data.table(
  parent = c("a", "b", "b", "b", "c", "g"),
  child = c("b", "c", "d", "e", "f", "h"),
  type = c("part_of", "part_of", "part_of", "is_a", "is_a", "part_of")
)

weights_v1 <- c("part_of" = 0.8, "is_a" = 0.6)

weights_v2 <- c("part_of" = 0.8, "is_a" = 0.8)

expected_wang_sim_v1 <- c(
  0.4698206,
  0,
  0.6247655,
  0.3995381,
  0.7960848,
  0,
  0.4278169,
  0,
  0.7641509,
  0.4767442,
  0.5901639,
  0,
  0.5575221,
  0,
  0,
  0,
  0.6428571,
  0,
  0.6428571,
  0.7641509,
  0,
  0.7422680,
  0.4767442,
  0,
  0.4805195,
  0,
  0.5575221,
  0
)

expected_wang_sim_v2 <- c(
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

expected_dt <- data.table::data.table(
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
  sim = expected_wang_sim_v1
)

expected_critval <- 0.7641509

## similarity matrices ---------------------------------------------------------

expect_error(
  current = calculate_wang_sim_mat(test_onto, weights = weights),
  info = "Wang ontology - correct error when column is missing"
)

# different values
results_v1 <- calculate_wang_sim_mat(test_onto_wang, weights = weights_v1)
critval_v1 <- calculate_critical_value(results_v1, alpha = 0.1)

expect_equivalent(
  current = rs_dense_to_upper_triangle(results_v1, 1L),
  target = expected_wang_sim_v1,
  info = "Wang similarity version 1 values",
  tolerance = 1e-6
)

expect_equivalent(
  current = critval_v1,
  target = expected_critval,
  info = "Wang similarity version 1 critical value",
  tolerance = 1e-6
)

# for the version 2 with different weights
results_v2 <- calculate_wang_sim_mat(test_onto_wang, weights = weights_v2)
critval_v2 <- calculate_critical_value(results_v2, alpha = 0.1)

expect_equivalent(
  current = rs_dense_to_upper_triangle(results_v2, 1L),
  target = expected_wang_sim_v2,
  info = "Wang similarity version 2 values",
  tolerance = 1e-6
)

expect_equivalent(
  current = critval_v2,
  target = expected_critval,
  info = "Wang similarity version 2 critical value",
  tolerance = 1e-6
)

## individual values -----------------------------------------------------------

expected_individual_res <- data.table::data.table(
  term1 = c("a", "a", "a", "b", "b", "e"),
  term2 = c("b", "e", "g", "e", "g", "g"),
  sims = c(0.6428571, 0.4805195, 0, 0.7422680, 0, 0)
)

individual_results <- calculate_wang_sim(
  terms = c("a", "b", "e", "g"),
  parent_child_dt = test_onto_wang,
  weights = weights_v1
)

expect_equal(
  current = individual_results,
  target = individual_results,
  tolerance = 1e-6,
  info = "Wang similarity - individual terms"
)

expect_error(
  current = calculate_wang_sim(
    terms = c("a", "b", "e", "g", "x"),
    parent_child_dt = test_onto_wang,
    weights = weights_v1
  ),
  info = "Wang similarity - individual terms error: wrong term"
)

expect_error(
  current = calculate_wang_sim(
    terms = c("a", "b", "e", "g"),
    parent_child_dt = test_onto,
    weights = weights_v1
  ),
  info = "Wang similarity - individual terms error: wrong parent_child_dt"
)

## class -----------------------------------------------------------------------

test_class <- ontology(test_onto, .verbose = FALSE)

expect_error(
  current = calculate_wang_sim_onto(test_class, .verbose = FALSE),
  info = "Wang ontology class - correct error when column is missing"
)

test_class <- ontology(test_onto_wang, .verbose = FALSE)

test_class <- calculate_wang_sim_onto(
  object = test_class,
  weights = weights_v1,
  .verbose = FALSE
)

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
  target = expected_wang_sim_v1,
  info = "Ontology class wang similarity - matrix version",
  tolerance = 1e-6
)

expect_equal(
  current = dt_res,
  target = expected_dt,
  info = paste(
    "Ontology class wang similarity - data.table version"
  ),
  tolerance = 1e-6
)

expect_equivalent(
  current = critval_class,
  target = expected_critval,
  info = "Ontology class wang similarity - critical value",
  tolerance = 1e-6
)
