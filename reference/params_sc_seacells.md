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
  knn = list()
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

- knn:

  List. Optional overrides for kNN parameters. See
  [`params_knn_defaults()`](params_knn_defaults.md) for available
  parameters: `k`, `knn_method`, `ann_dist`, `search_budget`, `n_trees`,
  `nn_max_iter`, `rho`, `delta`. Note: This function uses `k = 25L` as
  default.

## Value

A list with the SEACells parameters.
