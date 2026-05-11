# Wrapper function for parameters for SuperCell generation

Wrapper function for parameters for SuperCell generation

## Usage

``` r
params_sc_supercell(
  walk_length = 3L,
  graining_factor = 20,
  use_kernel = TRUE,
  k_ith = NULL,
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

- use_kernel:

  Boolean. Shall a kernel function akin to MAGIC be applied akin to the
  approach in SuperCell2, see Hérault, et al., bioRxiv, 2026 and van
  Dijk, et al., Cell, 2018.

- knn:

  List. Optional overrides for kNN parameters. See
  [`params_knn_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_defaults.md)
  for available parameters: `k`, `knn_method`, `ann_dist`,
  `search_budget`, `n_trees`, `delta`, `diversify_prob`, `ef_budget`,
  `m`, `ef_construction`, `ef_search`, `n_list` and `n_probe`.

- k_ith_neighbour:

  Optional integer. The k-ith neighbour to use for the kernel. Defaults
  to `k %/% 2`.

## Value

A list with the SuperCell parameters.
