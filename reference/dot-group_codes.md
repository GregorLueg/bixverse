# Dense 0-based codes for a grouping vector handed to Rust

Levels are sorted by code point (`method = "radix"`), so the codes do
not depend on the locale and match the Python front end. NA labels
error: cast to an integer they would reach Rust as `NA_integer_`.

## Usage

``` r
.group_codes(x, what)
```

## Arguments

- x:

  Atomic vector or list column of group labels.

- what:

  String. Name of the grouping for the error message.

## Value

Integer vector of 0-based codes, one per element of `x`.
