# Wrapper function for parameters for HotSpot

`weighted_graph` controls how the kNN distances become edge weights. The
default of `FALSE` follows the reference implementation: the distances
only decide who is a neighbour and every retained edge weighs one. Set
it to `TRUE` for the Gaussian kernel, whose width is the
`ceil(k / neighborhood_factor)`-th neighbour distance.

## Usage

``` r
params_sc_hotspot(
  model = c("danb", "normal", "bernoulli"),
  normalise = TRUE,
  weighted_graph = FALSE,
  neighborhood_factor = 3,
  knn = list()
)
```

## Arguments

- model:

  String. Model to use for modelling the GEX. One of
  `c("danb", "bernoulli", "normal")`. Defaults to `"danb"`.

- normalise:

  Boolean. Shall the data be normalised. Defaults to `TRUE`.

- weighted_graph:

  Boolean. Shall the Gaussian kernel be applied to the neighbour
  distances. Defaults to `FALSE`.

- neighborhood_factor:

  Float. Kernel width is the `ceil(k / neighborhood_factor)`-th
  neighbour distance. Only read when `weighted_graph = TRUE`. Defaults
  to `3`.

- knn:

  List. Optional overrides for kNN parameters. See
  [`params_knn_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_defaults.md)
  for available parameters: `k`, `knn_method`, `ann_dist`,
  `search_budget`, `n_trees`, `delta`, `diversify_prob`, `ef_budget`,
  `extract_knn`, `m`, `ef_construction`, `ef_search`, `n_list` and
  `n_probe`.

## Value

A list with the HotSpot parameters.

## References

DeTomaso and Yosef, Cell Systems, 2021
