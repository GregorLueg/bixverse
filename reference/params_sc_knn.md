# Parameters for single cell kNN searches

Parameters for single cell kNN searches

## Usage

``` r
params_sc_knn(
  k = 15L,
  knn_method = c("kmknn", "hnsw", "annoy", "nndescent", "ivf", "exhaustive"),
  ann_dist = c("euclidean", "cosine"),
  n_trees = 50L,
  search_budget = NULL,
  delta = 0.001,
  diversify_prob = 0,
  ef_budget = NULL,
  m = 16L,
  ef_construction = 200L,
  ef_search = 100L,
  n_list = NULL,
  n_probe = NULL
)
```

## Arguments

- k:

  Integer. Number of neighbours. Defaults to `15L`.

- knn_method:

  String. Which method to use for the approximate nearest neighbour
  search. `"nndescent"` hands back the graph the NN-Descent built, with
  no index and no query pass. Rows the descent left short are re-queried
  exactly. One of
  `c("kmknn", "hnsw", "annoy", "nndescent", "ivf", "exhaustive")`.
  Defaults to `"kmknn"`.

- ann_dist:

  String. Distance metric to use. One of `c("euclidean", "cosine")`.
  Defaults to `"euclidean"`.

- n_trees:

  Integer. Annoy param: number of trees. Defaults to `50L`.

- search_budget:

  Integer or `NULL`. Annoy param: optional search budget per tree. If
  `NULL`, defaults to `n_trees * k * 20L` internally. Defaults to
  `NULL`.

- delta:

  Numeric. NNDescent param: early termination criterion. Defaults to
  `0.001`.

- diversify_prob:

  Numeric. NNDescent param: diversification probability. Only used where
  a query runs against an NN-Descent index, not for the graph of a self
  kNN. Defaults to `0.0`.

- ef_budget:

  Integer or `NULL`. NNDescent param: optional query budget. Only used
  where a query runs against an NN-Descent index, not for the graph of a
  self kNN. Defaults to `NULL`.

- m:

  Integer. HNSW param: number of connections between layers. Defaults to
  `16L`.

- ef_construction:

  Integer. HNSW param: size of the dynamic candidate list during
  construction. Defaults to `200L`.

- ef_search:

  Integer. HNSW param: size of the candidate list at query time. Higher
  values improve recall at the cost of speed. Defaults to `100L`.

- n_list:

  Integer or `NULL`. IVF param: number of clusters to generate. If
  `NULL`, defaults to `sqrt(n)` internally. Defaults to `NULL`.

- n_probe:

  Integer or `NULL`. IVF param: number of clusters to query. If `NULL`,
  defaults to `sqrt(n_list)` internally. Defaults to `NULL`.

## Value

A named list with the following elements:

- k - Integer. Number of neighbours. Defaults to `15L`.

- knn_method - String. Which method to use for the approximate nearest
  neighbour search. `"nndescent"` hands back the graph the NN-Descent
  built, with no index and no query pass. Rows the descent left short
  are re-queried exactly. One of
  `c("kmknn", "hnsw", "annoy", "nndescent", "ivf", "exhaustive")`.
  Defaults to `"kmknn"`.

- ann_dist - String. Distance metric to use. One of
  `c("euclidean", "cosine")`. Defaults to `"euclidean"`.

- n_trees - Integer. Annoy param: number of trees. Defaults to `50L`.

- search_budget - Integer or `NULL`. Annoy param: optional search budget
  per tree. If `NULL`, defaults to `n_trees * k * 20L` internally.
  Defaults to `NULL`.

- delta - Numeric. NNDescent param: early termination criterion.
  Defaults to `0.001`.

- diversify_prob - Numeric. NNDescent param: diversification
  probability. Only used where a query runs against an NN-Descent index,
  not for the graph of a self kNN. Defaults to `0.0`.

- ef_budget - Integer or `NULL`. NNDescent param: optional query budget.
  Only used where a query runs against an NN-Descent index, not for the
  graph of a self kNN. Defaults to `NULL`.

- m - Integer. HNSW param: number of connections between layers.
  Defaults to `16L`.

- ef_construction - Integer. HNSW param: size of the dynamic candidate
  list during construction. Defaults to `200L`.

- ef_search - Integer. HNSW param: size of the candidate list at query
  time. Higher values improve recall at the cost of speed. Defaults to
  `100L`.

- n_list - Integer or `NULL`. IVF param: number of clusters to generate.
  If `NULL`, defaults to `sqrt(n)` internally. Defaults to `NULL`.

- n_probe - Integer or `NULL`. IVF param: number of clusters to query.
  If `NULL`, defaults to `sqrt(n_list)` internally. Defaults to `NULL`.
