# Wrapper function for parameters for neighbour identification in single cell

Wrapper function for parameters for neighbour identification in single
cell

## Usage

``` r
params_sc_neighbours(
  k = 15L,
  knn_algorithm = c("annoy", "hnsw", "nndescent"),
  ann_dist = c("cosine", "euclidean"),
  n_trees = 100L,
  search_budget = 100L,
  max_iter = 25L,
  rho = 1,
  delta = 0.001,
  full_snn = FALSE,
  pruning = 1/15,
  snn_similarity = c("rank", "jaccard")
)
```

## Arguments

- k:

  Integer. Number of neighbours to return.

- knn_algorithm:

  String. One of `c("annoy", "hnsw", "nndescent")`. Defaults to
  `"annoy"`.

- ann_dist:

  String. One of `c("cosine", "euclidean")`. The distance metric to be
  used for the approximate neighbour search.

- n_trees:

  Integer. Number of trees to use for the `annoy` algorithm.

- search_budget:

  Integer. Search budget per tree for the `annoy` algorithm.

- max_iter:

  Integer. Maximum iterations for the `nndescent` algorithm.

- rho:

  Numeric. Sampling rate for the `nndescent` algorithm.

- delta:

  Numeric. Early termination criterium for the `nndescent` algorithm.

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
  lower-ranked (closer) shared neighbors result in higher edge weights
  Both methods produce weights normalised to the range `[0, 1]`.

## Value

A list with the neighbour parameters.
