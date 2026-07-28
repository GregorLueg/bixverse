# Resolve the target matrix for a BulkCoExp method

Returns the pre-processed matrix if present, otherwise falls back to the
raw data with a warning. Used at the top of every bulk co-expression
method to remove duplicated preamble.

## Usage

``` r
.get_bulk_target_mat(object, .verbose = TRUE)
```

## Arguments

- object:

  The class, see
  [`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md).

- .verbose:

  Boolean. If `FALSE`, suppresses the fallback warning.

## Value

Numeric matrix.
