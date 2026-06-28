# Parameters for single cell kNN searches

Parameters for single cell kNN searches

## Usage

``` r
params_sc_knn(
  k = 15L,
  knn_method = "kmknn",
  ann_dist = "euclidean",
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
  search. Defaults to `"kmknn"`. One of
  `c("kmknn", "hnsw", "annoy", "nndescent", "ivf", "exhaustive")`.

- ann_dist:

  String. Distance metric to use. Defaults to `"euclidean"`. One of
  `c("cosine", "euclidean")`.

- n_trees:

  Integer. Annoy param: number of trees. Defaults to `50L`.

- search_budget:

  Integer or `NULL`. Annoy param: optional search budget per tree. If
  `NULL`, defaults to `n_trees * k * 20L` internally.

- delta:

  Numeric. NNDescent param: early termination criterion. Defaults to
  `0.001`.

- diversify_prob:

  Numeric. NNDescent param: diversification probability applied at the
  end of index construction. Defaults to `0.0`.

- ef_budget:

  Integer or `NULL`. NNDescent param: optional query budget. Higher
  values improve recall at the cost of speed.

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
  `NULL`, defaults to `sqrt(n)` internally.

- n_probe:

  Integer or `NULL`. IVF param: number of clusters to query. If `NULL`,
  defaults to `sqrt(n_list)` internally.

## Value

A list with the kNN parameters.
