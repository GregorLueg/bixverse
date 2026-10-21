# Run fast clustering on a SingleCells object

Coarsens the cells with k-means on the chosen embedding, builds a graph
over the centroids (a kNN graph with one of the symmetrisations of
[`params_knn_graph_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_graph_defaults.md),
or the cell-mediated graph) and runs Louvain or Leiden on it. The
community of each centroid is handed back to its cells. With
`n_seeds > 1`, every resolution is run once per seed, the labels come
from the seed with the lowest mean conductance and the stability across
the seeds is returned.

## Usage

``` r
fast_cluster_sc(
  object,
  embd_to_use = "pca",
  no_embd_to_use = NULL,
  resolutions = c(2, 1, 0.5),
  fc_params = params_sc_fast_cluster(),
  n_seeds = 1L,
  return_kmeans = FALSE,
  seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells` class.

- embd_to_use:

  String. Embedding name. Defaults to `"pca"`.

- no_embd_to_use:

  Optional integer. Number of dimensions to keep.

- resolutions:

  Numeric vector. Modularity resolutions.

- fc_params:

  List. Output of
  [`params_sc_fast_cluster()`](https://gregorlueg.github.io/bixverse/reference/params_sc_fast_cluster.md).

- n_seeds:

  Integer. Number of seeds per resolution. `1L` runs a single pass
  without stability statistics.

- return_kmeans:

  Boolean. Return k-means assignments and centroids.

- seed:

  Integer. Reproducibility.

- .verbose:

  Boolean or integer. Controls verbosity and returns run times. `FALSE`
  -\> quiet, `TRUE` or `1L` -\> normal verbosity, `2L` -\> detailed
  verbosity.

## Value

`SingleCellFastClusters` S3 object with:

- memberships:

  data.table with `cell_idx` and one column per resolution
  (`res_<value>`).

- stats:

  data.table with `resolution`, `mean_ari`, `median_ari`,
  `mean_conductance`, `median_conductance` and `mean_n_comms`, or `NULL`
  for `n_seeds = 1L`.

- k_means_cluster:

  Integer vector of k-means assignments, or `NULL`.

- centroids:

  Numeric matrix of centroids, or `NULL`.

- resolutions:

  Resolutions used.

with `cell_indices` stored as an attribute (0-indexed).

## References

Wei, et al., PLoS Comput Biol, 2022.

## Examples

``` r
# k-means centroids, Louvain on top, memberships back to the cells
sc <- demo_single_cells()
res <- fast_cluster_sc(
  sc,
  resolutions = c(1.0, 0.5),
  fc_params = params_sc_fast_cluster(n_centroids = 30L),
  .verbose = FALSE
)
head(get_data(res))
#>    cell_idx res_1 res_0.5
#>       <int> <int>   <int>
#> 1:        1     2       0
#> 2:        2     1       0
#> 3:        3     0       0
#> 4:        4     2       0
#> 5:        5     1       0
#> 6:        6     0       0

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
