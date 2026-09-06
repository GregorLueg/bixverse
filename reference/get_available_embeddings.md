# Get the available embeddings

Returns the available embedding as names from the class. This function
is used for the single cell-related classes and methods.

## Usage

``` r
get_available_embeddings(x, ...)

## S7 method for class <bixverse::MetaCells>
get_available_embeddings(x, ...)

# S3 method for class 'ScCache'
get_available_embeddings(x, ...)

## S7 method for class <bixverse::SingleCells>
get_available_embeddings(x, ...)

## S7 method for class <bixverse::SingleCellsSubset>
get_available_embeddings(x, ...)
```

## Arguments

- x:

  An object to get embedding from

- ...:

  Other parameters.

## Value

Get the names of the available embeddings.

## Examples

``` r
# what is in the cache to plot against
sc <- demo_single_cells()
get_available_embeddings(sc)
#> [1] "pca"

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
