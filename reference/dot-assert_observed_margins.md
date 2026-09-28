# Assert that no row or column is entirely missing

Checked on the R side so the error names 1-based indices.

## Usage

``` r
.assert_observed_margins(x)
```

## Arguments

- x:

  Numeric matrix.

## Value

Invisibly `x`, or an error naming the empty rows or columns.
