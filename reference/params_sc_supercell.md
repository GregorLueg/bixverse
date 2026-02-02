# Wrapper function for parameters for SuperCell generation

Wrapper function for parameters for SuperCell generation

## Usage

``` r
params_sc_supercell(
  walk_length = 3L,
  graining_factor = 20,
  linkage_dist = c("complete", "average"),
  knn = list()
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

- knn:

  List. Optional overrides for kNN parameters. See
  [`params_knn_defaults()`](params_knn_defaults.md) for available
  parameters: `k`, `knn_method`, `ann_dist`, `search_budget`, `n_trees`,
  `delta`, `diversify_prob`, `ef_budget`, `m`, `ef_construction`, and
  `ef_search`.

## Value

A list with the SuperCell parameters.
