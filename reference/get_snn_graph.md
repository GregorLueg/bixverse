# Get the sNN graph

Returns the shared nearest neighbour graph from the object. This
function is used for the single cell-related classes and methods.

## Usage

``` r
get_snn_graph(x, ...)

## S7 method for class <bixverse::MetaCells>
get_snn_graph(x, ...)

# S3 method for class 'ScCache'
get_snn_graph(x, ...)

## S7 method for class <bixverse::SingleCells>
get_snn_graph(x, ...)
```

## Arguments

- x:

  An object to get the sNN graph from.

- ...:

  Other parameters.

## Value

The igraph that has the shared nearest neighbours.
