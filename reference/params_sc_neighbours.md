# Wrapper function for parameters for neighbour identification in single cell

Wrapper function for parameters for neighbour identification in single
cell

## Usage

``` r
params_sc_neighbours(
  full_snn = TRUE,
  pruning = 0.0833333333333333,
  snn_similarity = c("jaccard", "rank"),
  knn = list()
)
```

## Arguments

- full_snn:

  Boolean. Shall the full shared nearest neighbour graph be generated
  that generates edges between all cells instead of between only
  neighbours. Defaults to `TRUE`.

- pruning:

  Numeric. Weights below this threshold will be set to 0 in the
  generation of the sNN graph. Seurat uses for example `1/15` with
  `k = 20`. As the default k is set to 15, we set it to `1/12`. Track
  this against `k` rather than leaving it: the threshold is a share of
  the neighbourhood, so the same value prunes far harder at a larger
  `k`. Over-pruning fails quietly, in that you still get a clustering,
  but cells left with too few shared neighbours drop out as singleton
  communities, which then show up downstream as one-cell clusters with
  inflated
  [`run_paga_sc()`](https://gregorlueg.github.io/bixverse/reference/run_paga_sc.md)
  connectivities. Defaults to `0.08333333333333333`.

- snn_similarity:

  String. The Jaccard similarity calculates the Jaccard between the
  neighbours, whereas the rank method calculates edge weights based on
  the ranking of shared neighbours. For the rank method, the weight is
  determined by finding the shared neighbour with the lowest combined
  rank across both cells, where lower-ranked (closer) shared neighbours
  result in higher edge weights Both methods produce weights normalised
  to the range `[0, 1]`. One of `c("jaccard", "rank")`. Defaults to
  `"jaccard"`.

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

- full_snn - Boolean. Shall the full shared nearest neighbour graph be
  generated that generates edges between all cells instead of between
  only neighbours. Defaults to `TRUE`.

- pruning - Numeric. Weights below this threshold will be set to 0 in
  the generation of the sNN graph. Seurat uses for example `1/15` with
  `k = 20`. As the default k is set to 15, we set it to `1/12`. Track
  this against `k` rather than leaving it: the threshold is a share of
  the neighbourhood, so the same value prunes far harder at a larger
  `k`. Over-pruning fails quietly, in that you still get a clustering,
  but cells left with too few shared neighbours drop out as singleton
  communities, which then show up downstream as one-cell clusters with
  inflated
  [`run_paga_sc()`](https://gregorlueg.github.io/bixverse/reference/run_paga_sc.md)
  connectivities. Defaults to `0.08333333333333333`.

- snn_similarity - String. The Jaccard similarity calculates the Jaccard
  between the neighbours, whereas the rank method calculates edge
  weights based on the ranking of shared neighbours. For the rank
  method, the weight is determined by finding the shared neighbour with
  the lowest combined rank across both cells, where lower-ranked
  (closer) shared neighbours result in higher edge weights Both methods
  produce weights normalised to the range `[0, 1]`. One of
  `c("jaccard", "rank")`. Defaults to `"jaccard"`.

- The elements of
  [`params_knn_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_knn_defaults.md),
  overridden by `knn`, spliced in at this position.
