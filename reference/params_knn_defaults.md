# Helper function to generate kNN defaults

This function generates various sensible default parameters for all of
the different approximate nearest neighbours that are available within
this package.

## Usage

``` r
params_knn_defaults()
```

## Value

A list with default parameters for kNN searches. Following parameters:

- k - Number of neighbours. Defaults to `15L`.

- knn_method - Which of method to use for the approximate nearest
  neighbour search. Defaults to `"hnsw"`. The implementations are:
  `c("hnsw", "annoy", "nndescent", "ivf", "exhaustive")`.

- ann_dist - Which distance metric to use for the approximate nearest
  neighbour search. Defaults to `"cosine"`. The implementations are
  `c("cosine", "euclidean")`.

- n_trees - Annoy param: number of trees to generate for Annoy. Defaults
  to `50L`.

- search_budget - Annoy param: optional search budget per tree for
  Annoy. If not provided, it will default to `n_tree * k * 20L`.

- diversify_prob - NNDescent param: diversification probability for the
  NNDescent index. This will diversify the index at the end and identify
  potentiall better edges. Defaults to `0.0`.

- delta - NNDescent param: early termination criterium for NNDescent.
  Defaults to `0.001`.

- ef_budget - NNDescent param: optional query budget parameter. Can
  accelerate querying, but at the cost of Recall.

- m - HNSW param: number of connections between layers for HNSW.
  Defaults to `16L`.

- ef_construction - HNSW param: size of dynamic candidate list during
  construction. Defaults to `200L`.

- ef_search - HNSW param: size of candidate list (higher = better
  recall, slower). Defaults to `100L`.

- n_list - IVF param: number of clusters/centroids to generate. Defaults
  to `NULL` (sqrt(n) n_lists will be generated in this case).

- n_probe - IVF param: number of clusters/centroids to query Defaults to
  `NULL` (sqrt(n_lists) clusters will be queried in this case).
