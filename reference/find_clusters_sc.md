# Graph-based clustering of cells on the sNN graph

This function will apply Leiden clustering on the sNN graph with the
given resolution and add a column to the obs table.

## Usage

``` r
find_clusters_sc(object, res = 1, name = "leiden_clustering")
```

## Arguments

- object:

  `single_cell_exp` class.

- res:

  Numeric. The resolution parameter for
  [`igraph::cluster_leiden()`](https://r.igraph.org/reference/cluster_leiden.html).

- name:

  String. The name to add to the obs table in the DuckDB.

## Value

The object with added clustering in the obs table.
