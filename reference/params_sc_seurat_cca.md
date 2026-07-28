# Wrapper function for the Seurat CCA parameters

Wrapper function for the Seurat CCA parameters

## Usage

``` r
params_sc_seurat_cca(
  num_cc = 30L,
  dims = 30L,
  k_anchor = 5L,
  k_filter = 200L,
  k_score = 30L,
  k_weight = 100L,
  n_top_features = 200L,
  l2_norm = TRUE,
  sd = 1,
  knn = list(),
  pca = params_sc_pca()
)
```

## Arguments

- num_cc:

  Integer. Number of canonical correlation dimensions to compute for the
  anchor space. Defaults to `30L`. The effective rank used is
  `max(num_cc, dims)`.

- dims:

  Integer. Number of dimensions used for the anchor kNN queries and the
  size of the returned embedding. Defaults to `30L`.

- k_anchor:

  Integer. Neighbourhood size for the mutual nearest neighbour anchor
  search. Defaults to `5L`.

- k_filter:

  Integer. Neighbourhood size for the gene-space anchor filter. Defaults
  to `200L`.

- k_score:

  Integer. Neighbourhood size for the shared-neighbour anchor scoring.
  Defaults to `30L`.

- k_weight:

  Integer. Neighbourhood size for the kernel weights applied during the
  correction. Defaults to `100L`.

- n_top_features:

  Integer. Number of top-loading genes used for the gene-space anchor
  filter. Defaults to `200L`.

- l2_norm:

  Boolean. Shall the canonical correlation embedding be L2-normalised
  per cell. Defaults to `TRUE`.

- sd:

  Numeric. Bandwidth divisor of the Gaussian kernel used for the anchor
  weights. Defaults to `1.0`.

- knn:

  List. Optional overrides for kNN parameters. See
  [`params_knn_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_defaults.md)
  for available parameters: `k`, `knn_method`, `ann_dist`,
  `search_budget`, `n_trees`, `delta`, `diversify_prob`, `ef_budget`,
  `m`, `ef_construction`, `ef_search`, `n_list` and `n_probe`. Note that
  `k` is unused here, the neighbourhood sizes come from `k_anchor`,
  `k_filter`, `k_score` and `k_weight`.

- pca:

  Named list. Parameters to feed through to the optional recalculation
  of the PCA, see
  [`params_sc_pca()`](https://gregorlueg.github.io/bixverse/reference/params_sc_pca.md).

## Value

A list with the Seurat CCA parameters.

## References

Stuart, et al., Cell, 2019
