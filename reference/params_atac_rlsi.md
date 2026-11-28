# Wrapper function for the rLSI integration parameters

Parameters of
[`run_rlsi_sc()`](https://gregorlueg.github.io/bixverse/reference/run_rlsi_sc.md),
Signac's `FindIntegrationAnchors(reduction = "rlsi")` followed by
`IntegrateEmbeddings()`.

## Usage

``` r
params_atac_rlsi(
  dims = 2:30,
  n_comp = 30L,
  k_anchor = 5L,
  k_score = 30L,
  k_weight = 100L,
  sd = 1,
  l2_norm = TRUE,
  solver = c("auto", "exact", "randomised", "covariance"),
  knn = list()
)
```

## Arguments

- dims:

  Integer vector. Per-batch LSI components the anchors use, 1-based. The
  default is Signac's `dims = 2:30`. Defaults to `2:30`.

- n_comp:

  Integer. Components of the per-batch LSIs; must reach `max(dims)`.
  Defaults to `30L`.

- k_anchor:

  Integer. Neighbours for the anchors. Defaults to `5L`.

- k_score:

  Integer. Neighbours for the anchor scores. Defaults to `30L`.

- k_weight:

  Integer. Anchors per cell in the correction weights. Defaults to
  `100L`.

- sd:

  Numeric. Bandwidth divisor of the weights. Defaults to `1.0`.

- l2_norm:

  Boolean. Scale and L2 normalise the anchor spaces. Defaults to `TRUE`.

- solver:

  String. SVD solver of the per-batch LSIs. One of
  `c("auto", "exact", "randomised", "covariance")`. Defaults to
  `"auto"`.

- knn:

  List. Optional overrides for kNN parameters. See
  [`params_knn_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_defaults.md)
  for available parameters. Only the backend is used; the neighbour
  counts are `k_anchor`, `k_score` and `k_weight`. See
  [`params_knn_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_defaults.md)
  for the available elements. Defaults to
  [`list()`](https://rdrr.io/r/base/list.html).

## Value

A named list with the following elements:

- dims - Integer vector. Per-batch LSI components the anchors use,
  1-based. The default is Signac's `dims = 2:30`. Defaults to `2:30`.

- n_comp - Integer. Components of the per-batch LSIs; must reach
  `max(dims)`. Defaults to `30L`.

- k_anchor - Integer. Neighbours for the anchors. Defaults to `5L`.

- k_score - Integer. Neighbours for the anchor scores. Defaults to
  `30L`.

- k_weight - Integer. Anchors per cell in the correction weights.
  Defaults to `100L`.

- sd - Numeric. Bandwidth divisor of the weights. Defaults to `1.0`.

- l2_norm - Boolean. Scale and L2 normalise the anchor spaces. Defaults
  to `TRUE`.

- solver - String. SVD solver of the per-batch LSIs. One of
  `c("auto", "exact", "randomised", "covariance")`. Defaults to
  `"auto"`.

- The elements of
  [`params_knn_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_defaults.md),
  overridden by `knn`, spliced in at this position.

## References

Stuart, et al., Nat Methods, 2021.
