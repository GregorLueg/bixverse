# Wrapper function for parameters for SuperCell generation

Wrapper function for parameters for SuperCell generation

## Usage

``` r
params_sc_supercell(
  walk_length = 3L,
  graining_factor = 20,
  linkage_dist = c("complete", "average"),
  k = 5L,
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

- walk_length:

  Integer. Walk length for the Walktrap algorithm. Defaults to `3L`.

- graining_factor:

  Numeric. Graining level of data (proportion of number of single cells
  in the initial dataset to the number of metacells in the final
  dataset). Defaults to `20.0`. (One meta cell per 20 cells.)

- linkage_dist:

  String. Which type of distance metric to use for the linkage. Defaults
  to `"average"`.

- k:

  Integer. Number of neighbours to return. Defaults to `5L`.

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

A list with the SuperCell parameters.
