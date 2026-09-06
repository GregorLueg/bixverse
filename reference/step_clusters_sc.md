# Pipeline step: graph-based clustering

Wraps
[`find_clusters_sc()`](https://gregorlueg.github.io/bixverse/reference/find_clusters_sc.md)
as an `ScStep`.

## Usage

``` r
step_clusters_sc(
  cluster_algorithm = c("leiden", "louvain"),
  res = 1,
  name = "leiden_clustering",
  modality = c("rna", "adt", "wnn"),
  seed = 42L
)
```

## Arguments

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

An `ScStep`.

## Examples

``` r
# Leiden over the graph the neighbours step wrote
step_neighbours_sc() %>>% step_clusters_sc(res = 0.5)
#> <ScPipeline> 2 steps
#>   1. neighbours  embd_to_use = "pca", no_embd_to_use = NULL, modality = c("rna", "adt"), neighbours_params = <list>, seed = 42L, .verbose = TRUE
#>   2. clusters    cluster_algorithm = c("leiden", "louvain"), res = 0.5, name = "leiden_clustering", modality = c("rna", "adt", "wnn"), seed = 42L
```
