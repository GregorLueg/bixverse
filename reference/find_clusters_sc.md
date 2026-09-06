# Graph-based clustering of cells on the sNN graph

This function will apply Leiden clustering on the sNN graph with the
given resolution and add a column to the obs table.

## Usage

``` r
find_clusters_sc(
  object,
  cluster_algorithm = c("leiden", "louvain"),
  res = 1,
  name = "leiden_clustering",
  modality = c("rna", "adt", "wnn"),
  seed = 42L
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

- modality:

  String. On which modality to run the UMAP. One of
  `c("rna", "adt", "wnn")`. The two latter options are only available
  for multi-modal versions with the added data.

- seed:

  Integer. For reproducibility.

## Value

The object with added clustering in the obs table.

## Examples

``` r
# Leiden on the cached sNN graph
sc <- demo_single_cells()
sc <- find_clusters_sc(sc, res = 1.0, name = "clusters")
table(get_sc_obs(sc)$clusters)
#> 
#>   0   1   2 
#> 169 166 165 

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
