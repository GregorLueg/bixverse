# Remove the neighbour graph

Remove the neighbour graph

## Usage

``` r
remove_connectivities(x, ...)

## S7 method for class <bixverse::MetaCells>
remove_connectivities(x, ...)

# S3 method for class 'ScCache'
remove_connectivities(x, ...)

## S7 method for class <bixverse::SingleCells>
remove_connectivities(x, ...)

## S7 method for class <bixverse::SingleCellsSubset>
remove_connectivities(x, ...)
```

## Arguments

- x:

  An object from which to remove the neighbour graph.

- ...:

  Other parameters.

## Examples

``` r
# drop the cached neighbour graph
sc <- demo_single_cells()
sc <- remove_connectivities(sc)
is.null(get_connectivities(sc))
#> [1] TRUE

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
