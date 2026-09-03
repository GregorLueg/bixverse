# Wrapper function for the Seurat rPCA parameters

Wrapper function for the Seurat rPCA parameters

## Usage

``` r
params_sc_seurat_rpca(
  dims = 30L,
  k_anchor = 5L,
  k_score = 30L,
  k_weight = 100L,
  l2_norm = TRUE,
  sd = 1,
  knn = list(),
  pca = params_sc_pca()
)
```

## Arguments

- dims:

  Integer. Number of dimensions used for the per-batch PCA projections,
  the anchor kNN queries and the size of the returned embedding.
  Defaults to `30L`.

- k_anchor:

  Integer. Neighbourhood size for the mutual nearest neighbour anchor
  search. Defaults to `5L`.

- k_score:

  Integer. Neighbourhood size for the shared-neighbour anchor scoring.
  Defaults to `30L`.

- k_weight:

  Integer. Neighbourhood size for the kernel weights applied during the
  correction. Defaults to `100L`.

- l2_norm:

  Boolean. Shall the projected embeddings be L2-normalised per cell.
  Defaults to `TRUE`.

- sd:

  Numeric. Bandwidth divisor of the Gaussian kernel used for the anchor
  weights. Defaults to `1.0`.

- knn:

  List. Optional overrides for kNN parameters. See
  [`params_knn_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_defaults.md)
  for available parameters: `k`, `knn_method`, `ann_dist`,
  `search_budget`, `n_trees`, `delta`, `diversify_prob`, `ef_budget`,
  `extract_knn`, `m`, `ef_construction`, `ef_search`, `n_list` and
  `n_probe`. Note that `k` is unused here, the neighbourhood sizes come
  from `k_anchor`, `k_score` and `k_weight`.

- pca:

  Named list. Parameters to feed through to the optional recalculation
  of the PCA, see
  [`params_sc_pca()`](https://gregorlueg.github.io/bixverse/reference/params_sc_pca.md).

## Value

A list with the Seurat rPCA parameters. Note that rPCA has no `num_cc`,
`k_filter` or `n_top_features`, the gene-space anchor filter is CCA-only
in Seurat.

## References

Stuart, et al., Cell, 2019
