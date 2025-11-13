# Replace values in a metadata column

This function will update the values in a given metadata column based on
what you are providing in terms of replacement.

## Usage

``` r
update_metadata_values(object, column, replacement, ...)
```

## Arguments

- object:

  The underlying object, either [`bixverse::bulk_coexp`](bulk_coexp.md)
  or [`bixverse::bulk_dge`](bulk_dge.md).

- column:

  Character vector. The columns for which to replace the values.

- replacement:

  Named character vector. The values with which to replace the data.

- ...:

  Additional arguments to parse to the functions.

## Value

Returns the object with the respective metadata updated.
