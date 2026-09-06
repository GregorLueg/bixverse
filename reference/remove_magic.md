# Remove the MAGIC imputed layer

Remove the MAGIC imputed layer

## Usage

``` r
remove_magic(x, ...)

# S3 method for class 'ScCache'
remove_magic(x, ...)

## S7 method for class <bixverse::SingleCells>
remove_magic(x, ...)

## S7 method for class <bixverse::SingleCellsSubset>
remove_magic(x, ...)
```

## Arguments

- x:

  An object from which to remove the imputed layer.

- ...:

  Other parameters.

## Value

The object with the imputed layer dropped.

## Examples

``` r
# drop the imputed layer again
sc <- demo_single_cells()
sc <- run_magic_sc(sc, features = get_gene_names(sc)[1:5], .verbose = FALSE)
sc <- remove_magic(sc)
is.null(get_magic(sc))
#> [1] TRUE

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
