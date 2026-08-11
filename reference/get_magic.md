# Get the MAGIC imputed layer

Returns the `ScMagic` layer written by
[`run_magic_sc()`](https://gregorlueg.github.io/bixverse/reference/run_magic_sc.md).
This function is used for the single cell-related classes and methods.

## Usage

``` r
get_magic(x, ...)

# S3 method for class 'ScCache'
get_magic(x, ...)

## S7 method for class <bixverse::SingleCells>
get_magic(x, ...)

## S7 method for class <bixverse::SingleCellsSubset>
get_magic(x, ...)
```

## Arguments

- x:

  An object to get the imputed layer from.

- ...:

  Other parameters.

## Value

The `ScMagic` object, or `NULL` when nothing was imputed.
