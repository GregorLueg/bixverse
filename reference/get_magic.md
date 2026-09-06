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

## Examples

``` r
# the imputed layer, only the genes MAGIC was asked for
sc <- demo_single_cells()
sc <- run_magic_sc(sc, features = get_gene_names(sc)[1:5], .verbose = FALSE)
dim(get_magic(sc)$data)
#> [1] 500   5

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
