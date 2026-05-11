# Fast single cell clustering parameters

Fast single cell clustering parameters

## Usage

``` r
params_sc_fast_cluster(
  kmeans_iters = 100L,
  batch_size = 4096L,
  drift_threshold = 1e-04,
  lr_alpha = 1,
  full_snn = FALSE,
  pruning = NULL,
  snn_similarity = c("jaccard", "rank"),
  louvain_iters = 10L,
  knn = list(k = 5L)
)
```

## Arguments

- kmeans_iters:

  Integer. Number of iterations for k-means clustering.

- batch_size:

  Integer. Batch size for mini batch k-means clustering.

- drift_threshold:

  Numeric. The drift for the mini batch k-means clustering. If the
  centroid drift is below this, the mini batch k-means terminates.

- lr_alpha:

  Numeric. Learning rate alpha parameter for mini batch k-means.

- full_snn:

  Boolean. Shall the full shared nearest neighbour graph be generated
  that generates edges between all cells instead of between only
  neighbours.

- pruning:

  Optional numeric. Weights below this threshold will be set to 0 in the
  generation of the sNN graph. If not provided, defaults to
  `1 / ceil(k * 0.8)`.

- snn_similarity:

  String. One of `c("rank", "jaccard")`. The Jaccard similarity
  calculates the Jaccard between the neighbours, whereas the rank method
  calculates edge weights based on the ranking of shared neighbours. For
  the rank method, the weight is determined by finding the shared
  neighbour with the lowest combined rank across both cells, where
  lower-ranked (closer) shared neighbours result in higher edge weights
  Both methods produce weights normalised to the range `[0, 1]`.

- louvain_iters:

  Integer. Number of iterations for Louvain clustering.

- knn:

  List. Optional overrides for kNN parameters. See
  [`params_knn_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_defaults.md)
  for available parameters: `k`, `knn_method`, `ann_dist`,
  `search_budget`, `n_trees`, `delta`, `diversify_prob`, `ef_budget`,
  `m`, `ef_construction`, `ef_search`, `n_list` and `n_probe`. Sets the
  default `k = 5L`.

## Value

A named list with the single cell fast clustering parameters.
