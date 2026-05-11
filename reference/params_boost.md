# Wrapper function for Boost parameters

Wrapper function for Boost parameters

## Usage

``` r
params_boost(
  boost_rate = 0.25,
  replace = FALSE,
  resolution = 1,
  n_iters = 20L,
  p_thresh = 1e-07,
  voter_thresh = 0.9,
  fast_cluster = FALSE,
  normalisation = list(),
  hvg = list(),
  pca = list(),
  knn = list(k = 0L),
  fast_cluster_params = list()
)
```

## Arguments

- boost_rate:

  Numeric. Boosting rate for the algorithm. Must be between 0 and 1.
  Defaults to `0.25`.

- replace:

  Boolean. Whether to use replacement during boosting. Defaults to
  `FALSE`.

- resolution:

  Numeric. Resolution parameter for graph-based clustering. Higher
  values lead to more clusters. Defaults to `1.0`.

- n_iters:

  Integer. Number of iterations to run the algorithm. Defaults to `20L`.

- p_thresh:

  Numeric. P-value threshold for significance testing. Defaults to
  `1e-7`.

- voter_thresh:

  Numeric. Voter threshold across iterations. Proportion of iterations a
  cell must be assigned to a cluster to be considered a member. Must be
  between 0 and 1. Defaults to `0.9`.

- fast_cluster:

  Boolean. Shall fast Louvain clustering be applied, i.e., k-means
  clustering and use the centroids for kNN graph generation and Louvain
  clustering with then backpropagating the membership based on centroid
  proximity.

- normalisation:

  List. Optional overrides for normalisation parameters. See
  [`params_norm_doublets_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_norm_doublets_defaults.md)
  for available parameters: `log_transform`, `mean_center`,
  `normalise_variance`, `target_size`. Note: Boost uses different
  defaults (`log_transform = FALSE`, `mean_center = TRUE`,
  `normalise_variance = TRUE`, `target_size = NULL`).

- hvg:

  List. Optional overrides for highly variable gene selection
  parameters. See
  [`params_hvg_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_hvg_defaults.md)
  for available parameters: `min_gene_var_pctl`, `hvg_method`,
  `loess_span`, `clip_max`.

- pca:

  List. Optional overrides for PCA parameters. See
  [`params_pca_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_pca_defaults.md)
  for available parameters: `no_pcs`, `random_svd`.

- knn:

  List. Optional overrides for kNN parameters. See
  [`params_knn_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_defaults.md)
  for available parameters: `k`, `knn_method`, `ann_dist`,
  `search_budget`, `n_trees`, `delta`, `diversify_prob`, `ef_budget`,
  `m`, `ef_construction`, `ef_search`, `n_list` and `n_probe`. Note:
  this function defaults to `k = 0L` (automatic neighbour detection).

- fast_cluster_params:

  List. Optional overrides for the fast clustering parameters. Only
  relevant if `fast_cluster = TRUE`. See
  [`params_fast_cluster_default()`](https://gregorlueg.github.io/bixverse/reference/params_fast_cluster_default.md)
  for available parameters: `km_type`, `n_centroids`, `kmeans_iters` and
  `batch_size`.

## Value

A named list with all Boost parameters, combining defaults with any
user-specified overrides.

A list with the Boost parameters.
