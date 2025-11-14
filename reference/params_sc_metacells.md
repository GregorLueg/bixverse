# Wrapper function for parameters for meta cell generation

Wrapper function for parameters for meta cell generation

## Usage

``` r
params_sc_metacells(
  max_shared = 15L,
  target_no_metacells = 1000L,
  max_iter = 5000L,
  knn = list()
)
```

## Arguments

- max_shared:

  Integer. Maximum number of allowed shared neighbours for the meta cell
  to be considered. Defaults to `15L`.

- target_no_metacells:

  Integer. Target number of meta-cells to generate. Defaults to `1000L`.

- max_iter:

  Integer. Maximum number of iterations for the algorithm. Defaults to
  `5000L`.

- knn:

  List. Optional overrides for kNN parameters. See
  [`params_knn_defaults()`](params_knn_defaults.md) for available
  parameters: `k`, `knn_method`, `ann_dist`, `search_budget`, `n_trees`,
  `nn_max_iter`, `rho`, `delta`. Note: This function uses `k = 25L` and
  `ann_dist = "cosine"` as defaults.

## Value

A list with the metacell parameters.
