# Graph-based clustering of cells on the sNN graph

This function will apply Leiden clustering on the sNN graph with the
given resolution and add a column to the obs table.

## Usage

``` r
find_clusters_sc(
  object,
  cluster_algorithm = c("leiden", "louvain"),
  res = 1,
  name = "leiden_clustering"
)
```

## Arguments

- object:

  `SingleCells`, `MetaCells` (or potentially other) class.

- cluster_algorithm:

  String. One of `c("leiden", "louvain")`.

- res:

  Numeric. The resolution parameter for
  [`igraph::cluster_leiden()`](https://r.igraph.org/reference/cluster_leiden.html)
  or
  [`igraph::cluster_louvain()`](https://r.igraph.org/reference/cluster_louvain.html).

- name:

  String. The name to add to the obs table in the DuckDB.

## Value

The object with added clustering in the obs table.
