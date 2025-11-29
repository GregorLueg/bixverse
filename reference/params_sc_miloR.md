# Wrapper function for parameters for MiloR

Wrapper function for parameters for MiloR

## Usage

``` r
params_sc_miloR(
  prop = 0.2,
  k_refine = 20L,
  refinement_strategy = c("index", "approximate", "bruteforce"),
  index_type = c("annoy", "hnsw"),
  knn = list()
)
```

## Arguments

- prop:

  Numeric. Proportion of cells to sample as neighbourhood indices.
  Defaults to `0.2`. Must be in (0,1).

- k_refine:

  Integer. Number of neighbours to use for refinement. Defaults to
  `20L`.

- refinement_strategy:

  String. Strategy for refining sampled indices. One of
  `c("approximate", "bruteforce", "index")`. Defaults to `"index"`.

- index_type:

  String. Type of kNN index to use. One of `c("annoy", "hnsw")`.
  Defaults to `"annoy"`.

- knn:

  List. Optional overrides for kNN parameters. See
  [`params_knn_defaults()`](params_knn_defaults.md) for available
  parameters: `k`, `knn_method`, `ann_dist`, `search_budget`, `n_trees`,
  `nn_max_iter`, `rho`, `delta`. Note: `knn_method` cannot be
  `"nndescent"` for MiloR as it doesn't generate an index!

## Value

A list with the MiloR parameters.
