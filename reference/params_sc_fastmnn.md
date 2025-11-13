# Wrapper function for the fastMNN parameters

Wrapper function for the fastMNN parameters

## Usage

``` r
params_sc_fastmnn(
  k = 20L,
  sigma = 0.1,
  knn_method = c("annoy", "hnsw"),
  dist_metric = c("cosine", "euclidean"),
  annoy_n_trees = 100L,
  annoy_search_budget = 100L,
  cos_norm = TRUE,
  var_adj = TRUE,
  no_pcs = 30L,
  random_svd = TRUE
)
```

## Arguments

- k:

  Integer. Number of mutual nearest neighbours to identify. Defaults to
  `20L`.

- sigma:

  Numeric. Bandwidth of the Gaussian smoothing kernel (as proportion of
  space radius). Defaults to `0.1`.

- knn_method:

  String. One of `c("annoy", "hnsw")`. Defaults to `"annoy"`.

- dist_metric:

  String. One of `c("cosine", "euclidean")`. Defaults to `"cosine"`.

- annoy_n_trees:

  Integer. Number of trees for Annoy index. Defaults to `100L`.

- annoy_search_budget:

  Integer. Search budget per tree for Annoy. Defaults to `100L`.

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

## Value

A list with the fastMNN parameters.
