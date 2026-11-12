# Graph-based clustering of cells on the neighbour graph

This function will apply Leiden or Louvain on the neighbour graph with
the given resolution and add a column to the obs table. Runs in Rust via
[`leiden_clustering()`](https://gregorlueg.github.io/bixverse/reference/leiden_clustering.md);
the labels depend on `seed` only, `0` is the largest cluster.

## Usage

``` r
find_clusters_sc(
  object,
  cluster_algorithm = c("leiden", "louvain"),
  res = 1,
  name = NULL,
  modality = c("rna", "adt", "wnn"),
  seed = 42L,
  .verbose = TRUE,
  ...
)
```

## Arguments

- object:

  `SingleCells`, `MetaCells` (or potentially other) class.

- cluster_algorithm:

  String. One of `c("leiden", "louvain")`.

- res:

  Numeric. The modularity resolution parameter, see
  [`leiden_clustering()`](https://gregorlueg.github.io/bixverse/reference/leiden_clustering.md).

- name:

  Optional string. The name of the obs column to add. `NULL` gives
  `"leiden_clustering"` or `"louvain_clustering"`, following
  `cluster_algorithm`.

- modality:

  String. On which modality to cluster. One of `c("rna", "adt", "wnn")`.
  The two latter options are only available for multi-modal versions
  with the added data.

- seed:

  Integer. For reproducibility.

- .verbose:

  Boolean or integer. Controls verbosity.

- ...:

  Further parameters for
  [`leiden_clustering()`](https://gregorlueg.github.io/bixverse/reference/leiden_clustering.md):
  `n_iterations`, `beta`, `min_gain` and `local_moving`.

## Value

The object with added clustering in the obs table.

## Examples

``` r
# Leiden on the cached neighbour graph
sc <- demo_single_cells()
sc <- find_clusters_sc(sc, res = 1.0, name = "clusters")
#> Running leiden clustering at resolution 1.
table(get_sc_obs(sc)$clusters)
#> 
#>   0   1   2 
#> 169 166 165 

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
