# Default parameters for Harmony v2 batch correction

Default parameters for Harmony v2 batch correction

## Usage

``` r
params_sc_harmony_v2(
  k = NULL,
  sigma = 0.1,
  theta = 2,
  lambda = 1,
  block_size = 0.2,
  max_iter_kmeans = 4L,
  max_iter_harmony = 10L,
  epsilon_kmeans = 0.001,
  epsilon_harmony = 0.01,
  window_size = 3L,
  alpha = 0.2,
  tau = 0,
  batch_proportion_cutoff = 1e-05,
  use_dynamic_lambda = FALSE,
  kmeans = list()
)
```

## Arguments

- k:

  Integer or `NULL`. Number of clusters for k-means clustering. If not
  provided, it will be automatically determined as
  `min(round(N / 30), 100)`. Defaults to `NULL`.

- sigma:

  Numeric vector. Per-cluster diversity weights. Either a single value
  (broadcast to all clusters) or a vector of length k. Defaults to
  `0.1`.

- theta:

  Numeric vector. Per-variable diversity penalties. Either a single
  value (broadcast to all variables) or a vector of length equal to the
  number of batch variables. Defaults to `2.0`.

- lambda:

  Numeric vector. Ridge regression penalty for the linear model.
  Typically a single value that is broadcast to all design matrix
  columns. Ignored when `use_dynamic_lambda = TRUE`. Defaults to `1.0`.

- block_size:

  Numeric. Fraction of cells to update per block during optimisation
  (0.0-1.0). Lower values reduce memory usage but increase computation
  time. Defaults to `0.2`.

- max_iter_kmeans:

  Integer. Maximum number of k-means iterations per Harmony round.
  Defaults to `4L`.

- max_iter_harmony:

  Integer. Maximum number of Harmony outer iterations. Defaults to
  `10L`.

- epsilon_kmeans:

  Numeric. Convergence threshold for k-means clustering. Stops when the
  relative change in cluster assignments falls below this value.
  Defaults to `0.001`.

- epsilon_harmony:

  Numeric. Convergence threshold for Harmony. Stops when the relative
  change in the objective function falls below this value. Defaults to
  `0.01`.

- window_size:

  Integer. Number of previous iterations to consider when checking
  convergence. Defaults to `3L`.

- alpha:

  Numeric. Scaling factor for dynamic lambda estimation. Must be in (0,
  1). Only relevant when `use_dynamic_lambda = TRUE`. Defaults to `0.2`.

- tau:

  Numeric. Scaling factor for theta based on batch size. A value of 0
  disables batch-size scaling of theta. Defaults to `0.0`.

- batch_proportion_cutoff:

  Numeric. Cutoff for pruning batches with small proportions during
  ridge regression. Defaults to `1e-05`.

- use_dynamic_lambda:

  Boolean. If `TRUE`, lambda is estimated dynamically per cluster
  instead of using the fixed `lambda` value. Defaults to `FALSE`.

- kmeans:

  List. Optional overrides for the k-means clustering algorithm Possible
  parameters are `"k_means_iter"`, `"k_means_init"`, `"gemm"` and
  `"hamerly"`, see
  [`params_kmeans_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_kmeans_defaults.md).
  See
  [`params_kmeans_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_kmeans_defaults.md)
  for the available elements. Defaults to
  [`list()`](https://rdrr.io/r/base/list.html).

## Value

A named list with the following elements:

- k - Integer or `NULL`. Number of clusters for k-means clustering. If
  not provided, it will be automatically determined as
  `min(round(N / 30), 100)`. Defaults to `NULL`.

- sigma - Numeric vector. Per-cluster diversity weights. Either a single
  value (broadcast to all clusters) or a vector of length k. Defaults to
  `0.1`.

- theta - Numeric vector. Per-variable diversity penalties. Either a
  single value (broadcast to all variables) or a vector of length equal
  to the number of batch variables. Defaults to `2.0`.

- lambda - Numeric vector. Ridge regression penalty for the linear
  model. Typically a single value that is broadcast to all design matrix
  columns. Ignored when `use_dynamic_lambda = TRUE`. Defaults to `1.0`.

- block_size - Numeric. Fraction of cells to update per block during
  optimisation (0.0-1.0). Lower values reduce memory usage but increase
  computation time. Defaults to `0.2`.

- max_iter_kmeans - Integer. Maximum number of k-means iterations per
  Harmony round. Defaults to `4L`.

- max_iter_harmony - Integer. Maximum number of Harmony outer
  iterations. Defaults to `10L`.

- epsilon_kmeans - Numeric. Convergence threshold for k-means
  clustering. Stops when the relative change in cluster assignments
  falls below this value. Defaults to `0.001`.

- epsilon_harmony - Numeric. Convergence threshold for Harmony. Stops
  when the relative change in the objective function falls below this
  value. Defaults to `0.01`.

- window_size - Integer. Number of previous iterations to consider when
  checking convergence. Defaults to `3L`.

- alpha - Numeric. Scaling factor for dynamic lambda estimation. Must be
  in (0, 1). Only relevant when `use_dynamic_lambda = TRUE`. Defaults to
  `0.2`.

- tau - Numeric. Scaling factor for theta based on batch size. A value
  of 0 disables batch-size scaling of theta. Defaults to `0.0`.

- batch_proportion_cutoff - Numeric. Cutoff for pruning batches with
  small proportions during ridge regression. Defaults to `1e-05`.

- use_dynamic_lambda - Boolean. If `TRUE`, lambda is estimated
  dynamically per cluster instead of using the fixed `lambda` value.
  Defaults to `FALSE`.

- The elements of
  [`params_kmeans_defaults()`](https://gregorlueg.github.io/bixverse/reference/params_kmeans_defaults.md),
  overridden by `kmeans`, spliced in at this position.
