# Helper function to generate kNN defaults

Helper function to generate kNN defaults

## Usage

``` r
params_knn_defaults()
```

## Value

A list with default parameters for kNN searches. Following parameters:

- k - Number of neighbours. Defaults to `15L`.

- knn_method - Which of method to use for the approximate nearest
  neighbour search. Defaults to `"annoy"`. The implementations are:
  `c("annoy", "hnsw", "nndescent")`.

- ann_dist - Which distance metric to use for the approximate nearest
  neighbour search. Defaults to `"euclidean"`. The implementations are
  `c("euclidean", "cosine")`.

- search_budget - Search budget per tree for Annoy. Defaults to `100L`.

- n_trees - Number of trees to generate for Annoy. Defaults to `100L`.

- nn_max_iter - Maximum iterations for NNDescent. Defaults to `15L`.

- rho - Sampling rate for NNDescent. Defaults to `1.0`.

- delta - Early termination criterium for NNDescent. Defaults to
  `0.001`.
