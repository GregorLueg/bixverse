# Convert a graph into a SingleCellConnectivity

Moves a neighbour graph into Rust. An igraph is flattened into an edge
list (edge attribute `weight` if present, otherwise all weights 1) and
its graph attribute `method` is carried over. A `SingleCellConnectivity`
is returned as is.

## Usage

``` r
as_sc_connectivity(x, ...)

# S3 method for class 'igraph'
as_sc_connectivity(x, ...)

# S3 method for class 'SingleCellConnectivity'
as_sc_connectivity(x, ...)
```

## Arguments

- x:

  An undirected `igraph` or a `SingleCellConnectivity`.

- ...:

  Other parameters.

## Value

A `SingleCellConnectivity`.

## Examples

``` r
g <- igraph::make_graph("Zachary")
conn <- as_sc_connectivity(g)
conn
#> SingleCellConnectivity: 34 cells, 78 edges
#>   Method: NA
```
