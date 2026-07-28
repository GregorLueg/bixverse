# Apply testChoice rules by name

Validates the elements of a named list `x` against per-field
[`checkmate::testChoice()`](https://mllg.github.io/checkmate/reference/checkChoice.html)
choice sets. Fields whose names are not in `rules` are skipped. On
failure, returns an error string naming the first offending element and
appending an optional `hint`.

## Usage

``` r
apply_choice_rules(x, rules, label, hint = NULL)
```

## Arguments

- x:

  Named list of parameters to validate.

- rules:

  Named list mapping field name to the character vector of allowed
  choices.

- label:

  Short human-readable label used in the error message (e.g.
  `"MELD params"`).

- hint:

  Optional string appended to the error message. Defaults to `NULL` (no
  hint).

## Value

`TRUE` if all checked fields pass, otherwise a string of the form
`` "The element `<field>` in <label> is not one of the expected choices. <hint>" ``.
