# Run fast Louvain clustering on a SingleCells object

Runs k-means on the chosen embedding, builds a kNN graph on the
centroids, applies Louvain clustering and propagates memberships back to
the cells. Optionally runs a grid over multiple seeds and returns
stability statistics.

## Usage

``` r
fast_cluster_sc(
  object,
  embd_to_use = "pca",
  no_embd_to_use = NULL,
  resolutions = c(2, 1, 0.5),
  km_type = c("kmeans", "minibatch"),
  n_centroids = NULL,
  fc_params = params_sc_fast_cluster(),
  snn = TRUE,
  return_kmeans = FALSE,
  grid_search = FALSE,
  no_seeds = 10L,
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

  Numeric vector. Louvain resolutions.

- km_type:

  String. One of `c("kmeans", "minibatch")`. The former runs standard
  k-means, the latter a mini-batch version that can be useful for large
  data sets.

- n_centroids:

  Optional integer. Number of k-means centroids. Defaults to
  `sqrt(n_cells)` Rust-side if `NULL`.

- fc_params:

  List. Output of
  [`params_sc_fast_cluster()`](https://gregorlueg.github.io/bixverse/reference/params_sc_fast_cluster.md).

- snn:

  Boolean. Convert kNN to sNN.

- return_kmeans:

  Boolean. Return k-means assignments and centroids.

- grid_search:

  Boolean. Run multi-seed grid version.

- no_seeds:

  Integer. Number of additional seeds (only used when
  `grid_search = TRUE`).

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

  data.table of grid statistics, or `NULL`.

- k_means_cluster:

  Integer vector of k-means assignments, or `NULL`.

- centroids:

  Numeric matrix of centroids, or `NULL`.

- resolutions:

  Resolutions used.

with `cell_indices` stored as an attribute (0-indexed).
