# Wrapper function for parameters for MiloR

Wrapper function for parameters for MiloR

## Usage

``` r
params_sc_miloR(
  prop = 0.2,
  k_refine = 20L,
  refinement_strategy = c("index", "approximate", "bruteforce"),
  index_type = c("hnsw", "annoy", "nndescent", "ivf"),
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

  String. Type of kNN index to use. One of
  `c("hnsw", "annoy", "nndescent", "ivf")`. Defaults to `"hnsw"`.

- knn:

  List. Optional overrides for kNN parameters. See
  [`params_knn_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_defaults.md)
  for available parameters: `k`, `knn_method`, `ann_dist`,
  `search_budget`, `n_trees`, `delta`, `diversify_prob`, `ef_budget`,
  `m`, `ef_construction`, `ef_search`, `n_list` and `n_probe`. Note:
  `knn_method` cannot be `"exhaustive"` for MiloR as it doesn't generate
  an index!

## Value

A list with the MiloR parameters.
