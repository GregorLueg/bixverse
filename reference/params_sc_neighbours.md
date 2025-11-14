# Wrapper function for parameters for neighbour identification in single cell

Wrapper function for parameters for neighbour identification in single
cell

## Usage

``` r
params_sc_neighbours(
  full_snn = FALSE,
  pruning = 1/15,
  snn_similarity = c("rank", "jaccard"),
  knn = list()
)
```

## Arguments

- full_snn:

  Boolean. Shall the full shared nearest neighbour graph be generated
  that generates edges between all cells instead of between only
  neighbours.

- pruning:

  Numeric. Weights below this threshold will be set to 0 in the
  generation of the sNN graph.

- snn_similarity:

  String. One of `c("rank", "jaccard")`. The Jaccard similarity
  calculates the Jaccard between the neighbours, whereas the rank method
  calculates edge weights based on the ranking of shared neighbours. For
  the rank method, the weight is determined by finding the shared
  neighbour with the lowest combined rank across both cells, where
  lower-ranked (closer) shared neighbours result in higher edge weights
  Both methods produce weights normalised to the range `[0, 1]`.

- knn:

  List. Optional overrides for kNN parameters. See
  [`params_knn_defaults()`](params_knn_defaults.md) for available
  parameters: `k`, `knn_method`, `ann_dist`, `search_budget`, `n_trees`,
  `nn_max_iter`, `rho`, `delta`. Note: This function uses `k = 15L`,
  `ann_dist = "cosine"`, and `nn_max_iter = 25L` as defaults.

## Value

A list with the neighbour parameters.
