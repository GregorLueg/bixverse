# Replace values in a metadata column

This function will update the values in a given metadata column based on
what you are providing in terms of replacement.

## Usage

``` r
update_metadata_values(object, column, replacement, ...)
```

## Arguments

- object:

  The underlying object, either
  [`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md)
  or
  [`BulkDge()`](https://gregorlueg.github.io/bixverse/reference/BulkDge.md).

- column:

  Character vector. The columns for which to replace the values.

- replacement:

  Named character vector. The values with which to replace the data.

- ...:

  Additional arguments to parse to the functions.

## Value

Returns the object with the respective metadata updated.
