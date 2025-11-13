# Return the metadata

Getter function to extract the metadata from the
[`bulk_coexp()`](bulk_coexp.md) or [`bulk_dge()`](bulk_dge.md).

## Usage

``` r
get_metadata(object, ...)
```

## Arguments

- object:

  The underlying object, either [`bixverse::bulk_coexp`](bulk_coexp.md)
  or [`bixverse::bulk_dge`](bulk_dge.md).

- ...:

  Additional arguments to parse to the functions.

## Value

Returns the metadata stored in the class.
