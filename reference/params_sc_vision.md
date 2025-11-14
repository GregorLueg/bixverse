# Wrapper function for parameters for VISION with auto-correlation

Wrapper function for parameters for VISION with auto-correlation

## Usage

``` r
params_sc_vision(n_perm = 500L, n_cluster = 5L, knn = list())
```

## Arguments

- n_perm:

  Integer. Number of random gene sets to generate per cluster.

- n_cluster:

  Integer. Number of clusters for the random gene set clustering
  generation.

- knn:

  List. Optional overrides for kNN parameters. See
  [`params_knn_defaults()`](params_knn_defaults.md) for available
  parameters: `k`, `knn_method`, `ann_dist`, `search_budget`, `n_trees`,
  `nn_max_iter`, `rho`, `delta`. Note: This function uses `k = 15L` and
  `nn_max_iter = 25L` as defaults.

## Value

A list with the VISION parameters when you wish to use the
auto-correlation version.
