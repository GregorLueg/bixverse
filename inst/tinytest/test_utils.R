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

## tests -----------------------------------------------------------------------

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
