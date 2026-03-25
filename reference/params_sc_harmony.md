# Default parameters for Harmony batch correction

Default parameters for Harmony batch correction

## Usage

``` r
params_sc_harmony(
  k = NULL,
  sigma = 0.1,
  theta = 2,
  lambda = 1,
  block_size = 0.05,
  max_iter_kmeans = 20L,
  max_iter_harmony = 10L,
  epsilon_kmeans = 1e-05,
  epsilon_harmony = 1e-04,
  window_size = 2L
)
```

## Arguments

- k:

  Optional integer. Number of clusters for k-means clustering. If not
  provided, it will be automatically determined as
  `min(round(N / 30), 200)`.

- sigma:

  Numeric vector. Per-cluster diversity weights. Either a single value
  (broadcast to all clusters) or a vector of length k.

- theta:

  Numeric vector. Per-variable diversity penalties. Either a single
  value (broadcast to all variables) or a vector of length equal to the
  number of batch variables.

- lambda:

  Numeric vector. Ridge regression penalty for the linear model.
  Typically a single value that is broadcast to all design matrix
  columns.

- block_size:

  Numeric. Fraction of cells to update per block during optimisation
  (0.0-1.0). Lower values reduce memory usage but increase computation
  time.

- max_iter_kmeans:

  Integer. Maximum number of k-means iterations per Harmony round.

- max_iter_harmony:

  Integer. Maximum number of Harmony outer iterations.

- epsilon_kmeans:

  Numeric. Convergence threshold for k-means clustering. Stops when the
  relative change in cluster assignments falls below this value.

- epsilon_harmony:

  Numeric. Convergence threshold for Harmony. Stops when the relative
  change in the objective function falls below this value.

- window_size:

  Integer. Number of previous iterations to consider when checking
  convergence.

## Value

A list with the parameters.
