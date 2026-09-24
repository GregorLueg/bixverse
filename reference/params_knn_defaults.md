# Helper function to generate kNN defaults

This function generates various sensible default parameters for all of
the different approximate nearest neighbours that are available within
this package.

## Usage

``` r
params_knn_defaults()
```

## Value

A named list with the following elements:

- k - Integer. Number of neighbours. Defaults to `15L`.

- knn_method - String. Which method to use for the approximate nearest
  neighbour search. One of
  `c("kmknn", "hnsw", "annoy", "nndescent", "ivf", "exhaustive")`.
  Defaults to `"kmknn"`.

- ann_dist - String. Which distance metric to use for the approximate
  nearest neighbour search. One of `c("euclidean", "cosine")`. Defaults
  to `"euclidean"`.

- n_trees - Integer. Annoy param: number of trees to generate for Annoy.
  Defaults to `50L`.

- search_budget - Integer or `NULL`. Annoy param: optional search budget
  per tree for Annoy. If not provided, it will default to
  `n_tree * k * 20L`. Defaults to `NULL`.

- delta - Numeric. NNDescent param: early termination criterium for
  NNDescent. Defaults to `0.001`.

- diversify_prob - Numeric. NNDescent param: diversification probability
  for the NNDescent index. This will diversify the index at the end and
  identify potentially better edges. Defaults to `0.0`.

- ef_budget - Integer or `NULL`. NNDescent param: optional query budget
  parameter. Can accelerate querying, but at the cost of Recall.
  Defaults to `NULL`.

- extract_knn - Boolean. NNDescent param: hand back the graph the
  descent already built instead of beam searching it. Skips the query
  pass entirely, so it is much faster, at the cost of some recall. Rows
  the descent never filled come back padded with duplicate edges.
  Ignored by every other method. Defaults to `FALSE`.

- m - Integer. HNSW param: number of connections between layers for
  HNSW. Defaults to `16L`.

- ef_construction - Integer. HNSW param: size of dynamic candidate list
  during construction. Defaults to `200L`.

- ef_search - Integer. HNSW param: size of candidate list (higher =
  better recall, slower). Defaults to `100L`.

- n_list - Integer or `NULL`. IVF param: number of clusters/centroids to
  generate. `NULL` generates `sqrt(n)` lists. Defaults to `NULL`.

- n_probe - Integer or `NULL`. IVF param: number of clusters/centroids
  to query. `NULL` queries `sqrt(n_lists)` clusters. Defaults to `NULL`.
