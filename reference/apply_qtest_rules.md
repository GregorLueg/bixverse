# Apply qtest rules by name

Validates the elements of a named list `x` against per-field
[`checkmate::qtest()`](https://mllg.github.io/checkmate/reference/qassert.html)
patterns. Fields whose names are not in `rules` are skipped. On failure,
returns an error string naming the first offending element and appending
an optional `hint`.

## Usage

``` r
apply_qtest_rules(x, rules, label, hint = NULL)
```

## Arguments

- x:

  Named list of parameters to validate.

- rules:

  Named list mapping field name to a qtest pattern (or vector of
  patterns passed to `qtest`).

- label:

  Short human-readable label used in the error message (e.g.
  `"GSEA params"`).

- hint:

  Optional string appended to the error message to describe the expected
  types/ranges. Defaults to `NULL` (no hint).

## Value

`TRUE` if all checked fields pass, otherwise a string of the form
`` "The element `<field>` in <label> is invalid. <hint>" ``.
