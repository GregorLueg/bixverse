# Wrapper function for parameters for HotSpot

Wrapper function for parameters for HotSpot

## Usage

``` r
params_sc_hotspot(
  model = c("danb", "normal", "bernoulli"),
  normalise = TRUE,
  knn = list(ann_dist = "cosine")
)
```

## Arguments

- model:

  String. Model to use for modelling the GEX. One of
  `c("danb", "bernoulli", "normal")`. Defaults to `"danb"`.

- normalise:

  Boolean. Shall the data be normalised. Defaults to `TRUE`.

- knn:

  List. Optional overrides for kNN parameters. See
  [`params_knn_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_defaults.md)
  for available parameters: `k`, `knn_method`, `ann_dist`,
  `search_budget`, `n_trees`, `delta`, `diversify_prob`, `ef_budget`,
  `m`, `ef_construction`, `ef_search`, `n_list` and `n_probe`.

## Value

A list with the HotSpot parameters.
