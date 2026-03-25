# Wrapper function for the fastMNN parameters

Wrapper function for the fastMNN parameters

## Usage

``` r
params_sc_fastmnn(
  ndist = 3,
  cos_norm = TRUE,
  no_pcs = 30L,
  random_svd = TRUE,
  knn = list(k = 20L)
)
```

## Arguments

- ndist:

  Numeric. Number of median distances for the tricube kernel bandwidth.
  Defaults to `3.0`.

- cos_norm:

  Logical. Apply cosine normalisation before computing distances.
  Defaults to `TRUE`.

- no_pcs:

  Integer. Number of PCs to use for MNN calculations. Defaults to `30L`.

- random_svd:

  Logical. Use randomised SVD. Defaults to `TRUE`.

- knn:

  List. Optional overrides for kNN parameters. See
  [`params_knn_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_defaults.md)
  for available parameters: `k`, `knn_method`, `ann_dist`,
  `search_budget`, `n_trees`, `delta`, `diversify_prob`, `ef_budget`,
  `m`, `ef_construction`, `ef_search`, `n_list` and `n_probe`.

## Value

A list with the fastMNN parameters.
