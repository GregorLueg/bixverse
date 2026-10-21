# sNN graph accessors (deprecated)

**\[deprecated\]**

Renamed to
[`get_connectivities()`](https://gregorlueg.github.io/bixverse/reference/get_connectivities.md),
[`set_connectivities()`](https://gregorlueg.github.io/bixverse/reference/set_connectivities.md)
and
[`remove_connectivities()`](https://gregorlueg.github.io/bixverse/reference/remove_connectivities.md).
The cached graph can be the sNN, the UMAP fuzzy simplicial set or the
kNN union, see
[`params_knn_graph_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_graph_defaults.md).

## Usage

``` r
get_snn_graph(x, ...)

set_snn_graph(x, snn_graph, ...)

remove_snn_graph(x, ...)
```

## Arguments

- x:

  The object holding the graph.

- ...:

  Passed on to the replacement.

- snn_graph:

  Igraph. The graph to set.

## Value

See the replacements.
