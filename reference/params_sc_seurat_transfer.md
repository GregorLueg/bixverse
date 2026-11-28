# Wrapper function for the Seurat label transfer parameters

Parameters of
[`transfer_labels_seurat_sc()`](https://gregorlueg.github.io/bixverse/reference/transfer_labels_seurat_sc.md),
Seurat's `FindTransferAnchors(reduction = "cca")` followed by
`TransferData()`.

## Usage

``` r
params_sc_seurat_transfer(
  n_cc = 30L,
  k_anchor = 5L,
  k_score = 30L,
  scale_max = 10,
  l2_norm = TRUE,
  k_weight = 50L,
  sd = 1,
  knn = list()
)
```

## Arguments

- n_cc:

  Integer. Canonical correlation vectors. Defaults to `30L`.

- k_anchor:

  Integer. Neighbours for the anchors. Defaults to `5L`.

- k_score:

  Integer. Neighbours for the anchor scores. Defaults to `30L`.

- scale_max:

  Numeric. Clip of the scaled features. Defaults to `10.0`.

- l2_norm:

  Boolean. L2 normalise the CCA space. Defaults to `TRUE`.

- k_weight:

  Integer. Anchors per query cell in the weights. Defaults to `50L`.

- sd:

  Numeric. Bandwidth divisor of the weights. Defaults to `1.0`.

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

- n_cc - Integer. Canonical correlation vectors. Defaults to `30L`.

- k_anchor - Integer. Neighbours for the anchors. Defaults to `5L`.

- k_score - Integer. Neighbours for the anchor scores. Defaults to
  `30L`.

- scale_max - Numeric. Clip of the scaled features. Defaults to `10.0`.

- l2_norm - Boolean. L2 normalise the CCA space. Defaults to `TRUE`.

- k_weight - Integer. Anchors per query cell in the weights. Defaults to
  `50L`.

- sd - Numeric. Bandwidth divisor of the weights. Defaults to `1.0`.

- The elements of
  [`params_knn_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_defaults.md),
  overridden by `knn`, spliced in at this position.

## References

Stuart, et al., Cell, 2019.
