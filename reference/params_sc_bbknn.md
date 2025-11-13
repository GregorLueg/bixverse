# Wrapper function for the BBKNN parameters

Wrapper function for the BBKNN parameters

## Usage

``` r
params_sc_bbknn(
  neighbours_within_batch = 5L,
  knn_method = c("annoy", "hnsw"),
  ann_dist = c("cosine", "euclidean"),
  set_op_mix_ratio = 1,
  local_connectivity = 1,
  annoy_n_trees = 100L,
  search_budget = 100L,
  trim = NULL
)
```

## Arguments

- neighbours_within_batch:

  Integer. Number of neighbours to consider per batch.

- knn_method:

  String. One of `c("annoy", "hnsw")`. Defaults to `"annoy"`.

- ann_dist:

  String. One of `c("cosine", "euclidean")`. The distance metric to be
  used for the approximate neighbour search.

- set_op_mix_ratio:

  Numeric. Mixing ratio between union (1.0) and intersection (0.0).

- local_connectivity:

  Numeric. UMAP connectivity computation parameter, how many nearest
  neighbours of each cell are assumed to be fully connected.

- annoy_n_trees:

  Integer. Number of trees to use in the generation of the Annoy index.

- search_budget:

  Integer. Search budget per tree for the `annoy` algorithm.

- trim:

  Optional integer. Trim the neighbours of each cell to these many top
  connectivities. May help with population independence and improve the
  tidiness of clustering. If `NULL`, it defaults to
  `10 * neighbours_within_batch`.

## Value

A list with the BBKNN parameters.
