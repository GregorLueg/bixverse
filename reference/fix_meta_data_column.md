# Helper to fix meta-data columns to be R conform

This function will update the specified columns in the metadata of an
[`BulkDge()`](https://gregorlueg.github.io/bixverse/reference/BulkDge.md)
or
[`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md)
to be conform with R standard naming convetions. This is useful to do
before running DGE methods as they expect standardised names.

## Usage

``` r
fix_meta_data_column(object, col_names, ...)
```

## Arguments

- object:

  The underlying object, either
  [`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md)
  or
  [`BulkDge()`](https://gregorlueg.github.io/bixverse/reference/BulkDge.md).

- col_names:

  Character vector. The columns to fix.

- ...:

  Additional arguments to parse to the functions.

## Value

Returns the object with the respective metadata columns updated.
