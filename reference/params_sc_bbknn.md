# Wrapper function for the BBKNN parameters

Wrapper function for the BBKNN parameters

## Usage

``` r
params_sc_bbknn(
  neighbours_within_batch = 3L,
  set_op_mix_ratio = 1,
  local_connectivity = 1,
  trim = NULL,
  knn = list()
)
```

## Arguments

- neighbours_within_batch:

  Integer. Number of neighbours to consider per batch. Defaults to `3L`.

- set_op_mix_ratio:

  Numeric. Mixing ratio between union (1.0) and intersection (0.0).
  Defaults to `1.0`.

- local_connectivity:

  Numeric. UMAP connectivity computation parameter, how many nearest
  neighbours of each cell are assumed to be fully connected. Defaults to
  `1.0`.

- trim:

  Optional integer. Trim the neighbours of each cell to these many top
  connectivities. May help with population independence and improve the
  tidiness of clustering. If `NULL`, it defaults to
  `10 * neighbours_within_batch`.

- knn:

  List. Optional overrides for kNN parameters. See
  [`params_knn_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_defaults.md)
  for available parameters: `k`, `knn_method`, `ann_dist`,
  `search_budget`, `n_trees`, `delta`, `diversify_prob`, `ef_budget`,
  `m`, `ef_construction`, `ef_search`, `n_list` and `n_probe`.

## Value

A list with the BBKNN parameters.
