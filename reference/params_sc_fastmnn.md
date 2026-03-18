# Wrapper function for the fastMNN parameters

Wrapper function for the fastMNN parameters

## Usage

``` r
params_sc_fastmnn(
  sigma = 0.1,
  cos_norm = TRUE,
  var_adj = TRUE,
  no_pcs = 30L,
  random_svd = TRUE,
  knn = list(k = 20L)
)
```

## Arguments

- sigma:

  Numeric. Bandwidth of the Gaussian smoothing kernel (as proportion of
  space radius). Defaults to `0.1`.

- cos_norm:

  Logical. Apply cosine normalisation before computing distances.
  Defaults to `TRUE`.

- var_adj:

  Logical. Apply variance adjustment to avoid kissing effects. Defaults
  to `TRUE`.

- no_pcs:

  Integer. Number of PCs to use for MNN calculations. Defaults to `30L`.

- random_svd:

  Logical. Use randomised SVD. Defaults to `TRUE`.

- knn:

  List. Optional overrides for kNN parameters. See
  [`params_knn_defaults()`](params_knn_defaults.md) for available
  parameters: `k`, `knn_method`, `ann_dist`, `search_budget`, `n_trees`,
  `delta`, `diversify_prob`, `ef_budget`, `m`, `ef_construction`,
  `ef_search`, `n_list` and `n_probe`.

## Value

A list with the fastMNN parameters.
