# Wrapper function for the SEACells parameters

Wrapper function for the SEACells parameters

## Usage

``` r
params_sc_seacells(
  n_sea_cells,
  max_fw_iters = 50L,
  convergence_epsilon = 0.001,
  max_iter = 100L,
  min_iter = 10L,
  greedy_threshold = 20000L,
  graph_building = "union",
  pruning = FALSE,
  pruning_threshold = 1e-07,
  k = 25L,
  knn_method = c("annoy", "hnsw", "nndescent"),
  n_trees = 100L,
  search_budget = 100L,
  nn_max_iter = 15L,
  rho = 1,
  delta = 0.001
)
```

## Arguments

- n_sea_cells:

  Integer. Number of SEA cells to detect.

- max_fw_iters:

  Integer. Maximum iterations for the Franke-Wolfe algorithm. Defaults
  to `50L`.

- convergence_epsilon:

  Numeric. Convergence threshold. Algorithm stops when RSS change \<
  epsilon \* RSS(0). Defaults to `1e-3`.

- max_iter:

  Integer. Maximum iterations to run SEACells for. Defaults to `100L`.

- min_iter:

  Integer. Minimum iterations to run SEACells for. Defaults to `10L`.

- greedy_threshold:

  Integer. Maximum number of cells before defaulting to rapid random
  selection of archetypes. Defaults to `20000L`.

- graph_building:

  String. Graph building method. Defaults to `"union"`.

- pruning:

  Boolean. Shall tiny values be pruned during Franke-Wolfe updates. This
  will reduce memory pressure and can be a good option on large data
  sets. Defaults to `FALSE`.

- pruning_threshold:

  Float. If `pruning = TRUE` values below which threshold shall be
  pruned.

- k:

  Integer. Number of neighbours for the kNN algorithm. Defaults to
  `25L`.

- knn_method:

  String. One of `c("annoy", "hnsw", "nndescent")`. Defaults to
  `"annoy"`.

- n_trees:

  Integer. Number of trees for Annoy index. Defaults to `100L`.

- search_budget:

  Integer. Search budget during querying. Defaults to `100L`.

- nn_max_iter:

  Integer. Maximum iterations for NN Descent. Defaults to `15L`.

- rho:

  Numeric. Sampling rate for NN Descent. Defaults to `1.0`.

- delta:

  Numeric. Early termination criterion for NN Descent. Defaults to
  `0.001`.

## Value

A list with the SEACells parameters.
