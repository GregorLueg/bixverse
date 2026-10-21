# Set/add the neighbour graph

Set/add the neighbour graph

## Usage

``` r
set_connectivities(x, connectivities, ...)

## S7 method for class <bixverse::MetaCells>
set_connectivities(x, connectivities, ...)

# S3 method for class 'ScCache'
set_connectivities(x, connectivities, ...)

## S7 method for class <bixverse::SingleCells>
set_connectivities(x, connectivities, ...)

## S7 method for class <bixverse::SingleCellsSubset>
set_connectivities(x, connectivities, ...)
```

## Arguments

- x:

  An object to add the neighbour graph to.

- connectivities:

  Igraph. The neighbour graph for subsequent clustering, see
  [`find_neighbours_sc()`](https://gregorlueg.github.io/bixverse/reference/find_neighbours_sc.md).

- ...:

  Other parameters.

## Examples

``` r
# the neighbour graph taken out and put back
sc <- demo_single_cells()
g <- get_connectivities(sc)
sc <- set_connectivities(remove_connectivities(sc), g)
igraph::vcount(get_connectivities(sc))
#> [1] 500

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
