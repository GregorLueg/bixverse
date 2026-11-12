# Symmetrise a kNN into the neighbour graph

Symmetrise a kNN into the neighbour graph

## Usage

``` r
.knn_to_graph(knn, graph_params, .verbose = TRUE)
```

## Arguments

- knn:

  `SingleCellNearestNeighbour`. The kNN to symmetrise.

- graph_params:

  List. Carries the keys of
  [`params_knn_graph_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_graph_defaults.md).

- .verbose:

  Boolean or integer. Controls verbosity.

## Value

A `SingleCellConnectivity` holding the undirected, weighted Rust graph,
with the method in `method`.
