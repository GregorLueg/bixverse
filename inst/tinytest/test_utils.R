# util tests -------------------------------------------------------------------

## strings ---------------------------------------------------------------------

example_strings <- c(
  "Hello World",
  "camelCaseExample",
  "PascalCaseExample",
  "string-with-dashes",
  "string_with_underscores",
  "string.with.dots",
  "string@with#special$symbols%",
  " Multiple   Spaces",
  "Mixed-Case_and.Special@Characters!",
  "123numbers456",
  "ALLCAPS"
)

expected_strings <- c(
  "hello_world",
  "camel_case_example",
  "pascal_case_example",
  "string_with_dashes",
  "string_with_underscores",
  "string_with_dots",
  "string_with_special_symbols",
  "multiple_spaces",
  "mixed_case_and_special_characters",
  "123numbers456",
  "allcaps"
)

### tests ----------------------------------------------------------------------

new_strings <- to_snake_case(example_strings)

expect_equal(
  current = new_strings,
  target = expected_strings,
  info = "utils - to snake case transformer"
)

example_strings_2 <- c(example_strings, NA)

expect_error(
  current = to_snake_case(example_strings_2),
  info = "utils - snake case transformer: error with ignore_na = FALSE"
)

new_strings_2 <- to_snake_case(example_strings_2, ignore_na = TRUE)

expect_true(
  current = sum(is.na(new_strings_2)) == 1,
  info = "utils - snake case transformer: return NA with ignore_na = TRUE"
)

## knn transformations ---------------------------------------------------------

### rs_knn_mat_to_edge_list ----------------------------------------------------

# 3 samples, 2 neighbours, 0-indexed input
knn_mat <- matrix(
  as.integer(c(1, 2, 0, 2, 0, 1)),
  nrow = 3,
  ncol = 2,
  byrow = TRUE
)

result_one_idx <- rs_knn_mat_to_edge_list(knn_mat, one_index = TRUE)

expect_equal(
  current = length(result_one_idx),
  target = 3L * 2L * 2L,
  info = "edge list length is n_samples * k_neighbours * 2"
)
expect_equal(
  current = result_one_idx,
  target = as.integer(c(1, 2, 1, 3, 2, 1, 2, 3, 3, 1, 3, 2)),
  info = paste(
    "one_index = TRUE",
    "correctly converts 0-indexed matrix to 1-indexed edge list"
  )
)

result_zero_idx <- rs_knn_mat_to_edge_list(knn_mat, one_index = FALSE)

expect_equal(
  current = result_zero_idx,
  target = as.integer(c(0, 1, 0, 2, 1, 0, 1, 2, 2, 0, 2, 1)),
  info = "one_index = FALSE returns 0-indexed edge list unchanged"
)

### rs_knn_mat_to_edge_pairs ---------------------------------------------------

result_pairs <- rs_knn_mat_to_edge_pairs(knn_mat, one_index = TRUE)

expect_true(
  checkmate::test_list(result_pairs, names = "named"),
  info = "edge pairs returns a named list"
)
expect_equal(
  current = names(result_pairs),
  target = c("from", "to"),
  info = "edge pairs has correct names"
)
expect_equal(
  current = length(result_pairs$from),
  target = 3L * 2L,
  info = "from vector length is n_samples * k_neighbours"
)
expect_equal(
  current = length(result_pairs$to),
  target = 3L * 2L,
  info = "to vector length is n_samples * k_neighbours"
)
expect_equal(
  current = result_pairs$from,
  target = as.integer(c(1, 1, 2, 2, 3, 3)),
  info = "one_index = TRUE from indices are 1-based"
)
expect_equal(
  current = result_pairs$to,
  target = as.integer(c(2, 3, 1, 3, 1, 2)),
  info = "one_index = TRUE to indices are 1-based"
)

### consistency ----------------------------------------------------------------

# verify edge list and edge pairs are consistent
pairs_zero <- rs_knn_mat_to_edge_pairs(knn_mat, one_index = FALSE)
flat_zero <- rs_knn_mat_to_edge_list(knn_mat, one_index = FALSE)

expect_equal(
  current = as.integer(rbind(pairs_zero$from, pairs_zero$to)),
  target = flat_zero,
  info = "edge pairs and edge list produce consistent output"
)
