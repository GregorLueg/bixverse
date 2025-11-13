# Helper to fix meta-data columns to be R conform

This function will update the specified columns in the metadata of an
[`bixverse::bulk_dge`](bulk_dge.md) or
[`bixverse::bulk_coexp`](bulk_coexp.md) to be conform with R standard
naming convetions. This is useful to do before running DGE methods as
they expect standardised names.

## Usage

``` r
fix_meta_data_column(object, col_names, ...)
```

## Arguments

- object:

  The underlying object, either [`bixverse::bulk_coexp`](bulk_coexp.md)
  or [`bixverse::bulk_dge`](bulk_dge.md).

- col_names:

  Character vector. The columns to fix.

- ...:

  Additional arguments to parse to the functions.

## Value

Returns the object with the respective metadata columns updated.
