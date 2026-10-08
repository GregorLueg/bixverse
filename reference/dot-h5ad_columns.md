# Bring table columns into the types the h5ad writer encodes

The Rust side maps factors to categoricals, characters to string arrays,
and doubles, integers and logicals to plain arrays. This picks which of
those each column becomes:

- a character column with repeated values becomes a factor, which is
  what pandas would hold it as

- a column with nothing but missing values has no category to write and
  becomes the literal string `"NA"`

- integers and logicals with missing values are widened to double, the
  only one of the three with a representation for them (`NA_integer_` is
  `INT_MIN` on disk, a silently wrong number)

- anything else (dates, lists) is stringified rather than dropped

## Usage

``` r
.h5ad_columns(dt)
```

## Arguments

- dt:

  data.table. The columns to write, without the index.

## Value

A named list of factor, character, double, integer or logical vectors.
