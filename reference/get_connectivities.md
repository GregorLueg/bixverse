# Get the neighbour graph

Returns the neighbour graph the clustering runs on. Depending on the
`knn_sym` it was built with, see
[`params_knn_graph_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_graph_defaults.md),
this is the shared nearest neighbour graph, the UMAP fuzzy simplicial
set or the kNN union. The graph lives in Rust; convert it with
[`as_igraph()`](https://gregorlueg.github.io/bixverse/reference/as_igraph.md)
if you need an igraph.

## Usage

``` r
get_connectivities(x, ...)

## S7 method for class <bixverse::MetaCells>
get_connectivities(x, ...)

# S3 method for class 'ScCache'
get_connectivities(x, ...)

## S7 method for class <bixverse::SingleCells>
get_connectivities(x, ...)

## S7 method for class <bixverse::SingleCellsSubset>
get_connectivities(x, ...)
```

## Arguments

- x:

  An object to get the neighbour graph from.

- ...:

  Other parameters.

## Value

A `SingleCellConnectivity`, see
[`as_sc_connectivity()`](https://gregorlueg.github.io/bixverse/reference/as_sc_connectivity.md).

## Examples

``` r
# the neighbour graph the clustering methods run on
sc <- demo_single_cells()
get_connectivities(sc)
#> SingleCellConnectivity: 500 cells, 18888 edges
#>   Method: snn

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
