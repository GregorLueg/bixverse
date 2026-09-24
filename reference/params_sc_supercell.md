# Wrapper function for parameters for SuperCell generation

Wrapper function for parameters for SuperCell generation

## Usage

``` r
params_sc_supercell(
  walk_length = 3L,
  graining_factor = 20,
  use_kernel = TRUE,
  k_ith = NULL,
  max_support = NULL,
  knn = list()
)
```

## Arguments

- walk_length:

  Integer. Walk length for the Walktrap algorithm. Defaults to `3L`.

- graining_factor:

  Numeric. Graining level of data (proportion of number of single cells
  in the initial dataset to the number of metacells in the final
  dataset). (One meta cell per 20 cells.) Defaults to `20.0`.

- use_kernel:

  Boolean. Shall a kernel function akin to MAGIC be applied akin to the
  approach in SuperCell2, see Hérault, et al., bioRxiv, 2026 and van
  Dijk, et al., Cell, 2018. Defaults to `TRUE`.

- k_ith:

  Integer or `NULL`. The k-ith neighbour to use for the kernel. Defaults
  to `NULL`.

- max_support:

  Integer or `NULL`. Caps each cell's walk-probability vector to its top
  entries by mass, bounding memory at ~`max_support * n_cells` on large
  data. Makes the result an approximation. `NULL` (default) keeps the
  walks exact. Defaults to `NULL`.

- knn:

  List. Optional overrides for kNN parameters. See
  [`params_knn_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_defaults.md)
  for available parameters: `k`, `knn_method`, `ann_dist`,
  `search_budget`, `n_trees`, `delta`, `diversify_prob`, `ef_budget`,
  `extract_knn`, `m`, `ef_construction`, `ef_search`, `n_list` and
  `n_probe`. See
  [`params_knn_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_defaults.md)
  for the available elements. Defaults to
  [`list()`](https://rdrr.io/r/base/list.html).

## Value

A named list with the following elements:

- walk_length - Integer. Walk length for the Walktrap algorithm.
  Defaults to `3L`.

- graining_factor - Numeric. Graining level of data (proportion of
  number of single cells in the initial dataset to the number of
  metacells in the final dataset). (One meta cell per 20 cells.)
  Defaults to `20.0`.

- use_kernel - Boolean. Shall a kernel function akin to MAGIC be applied
  akin to the approach in SuperCell2, see Hérault, et al., bioRxiv, 2026
  and van Dijk, et al., Cell, 2018. Defaults to `TRUE`.

- k_ith - Integer or `NULL`. The k-ith neighbour to use for the kernel.
  Defaults to `NULL`.

- max_support - Integer or `NULL`. Caps each cell's walk-probability
  vector to its top entries by mass, bounding memory at
  ~`max_support * n_cells` on large data. Makes the result an
  approximation. `NULL` (default) keeps the walks exact. Defaults to
  `NULL`.

- The elements of
  [`params_knn_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_defaults.md),
  overridden by `knn`, spliced in at this position.
