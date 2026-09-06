# Set/add KNN

Set/add KNN

## Usage

``` r
set_snn_graph(x, snn_graph, ...)

## S7 method for class <bixverse::MetaCells>
set_snn_graph(x, snn_graph, ...)

# S3 method for class 'ScCache'
set_snn_graph(x, snn_graph, ...)

## S7 method for class <bixverse::SingleCells>
set_snn_graph(x, snn_graph, ...)

## S7 method for class <bixverse::SingleCellsSubset>
set_snn_graph(x, snn_graph, ...)
```

## Arguments

- x:

  An object to add the KNN data to.

- snn_graph:

  Igraph. The sNN graph for subsequent clustering.

- ...:

  Other parameters.

## Examples

``` r
# the sNN graph taken out and put back
sc <- demo_single_cells()
snn <- get_snn_graph(sc)
sc <- set_snn_graph(remove_snn_graph(sc), snn)
igraph::vcount(get_snn_graph(sc))
#> [1] 500

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
