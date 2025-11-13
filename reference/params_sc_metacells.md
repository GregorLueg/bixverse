# Wrapper function for parameters for meta cell generation

Wrapper function for parameters for meta cell generation

## Usage

``` r
params_sc_metacells(
  max_shared = 15L,
  target_no_metacells = 1000L,
  max_iter = 5000L,
  k = 25L,
  knn_method = c("annoy", "hnsw", "nndescent"),
  ann_dist = c("cosine", "euclidean"),
  n_trees = 100L,
  search_budget = 100L,
  nn_max_iter = 15L,
  rho = 1,
  delta = 0.001
)
```

## Arguments

- max_shared:

  Integer. Maximum number of allowed shared neighbours for the meta cell
  to be considered. Defaults to `15L`.

- target_no_metacells:

  Integer. Target number of meta-cells to generate. Defaults to `1000L`.

- max_iter:

  Integer. Maximum number of iterations for the algorithm. Defaults to
  `5000L`.

- k:

  Integer. Number of neighbours to return. Defaults to `25L`.

- knn_method:

  String. One of `c("annoy", "hnsw", "nndescent")`. Defaults to
  `"annoy"`.

- ann_dist:

  String. One of `c("cosine", "euclidean")`. The distance metric to be
  used for the approximate neighbour search. Defaults to `"cosine"`.

- n_trees:

  Integer. Number of trees to use for the `annoy` algorithm. Defaults to
  `100L`.

- search_budget:

  Integer. Search budget per tree for the `annoy` algorithm. Defaults to
  `100L`.

- nn_max_iter:

  Integer. Maximum iterations for NN Descent. Defaults to `15L`.

- rho:

  Numeric. Sampling rate for NN Descent. Defaults to `1.0`.

- delta:

  Numeric. Early termination criterion for NN Descent. Defaults to
  `0.001`.

## Value

A list with the metacell parameters.
