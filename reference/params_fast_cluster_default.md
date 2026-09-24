# Helper function to generate default parameters for the fast clustering for the doublet detection methods

Helper function to generate default parameters for the fast clustering
for the doublet detection methods

## Usage

``` r
params_fast_cluster_default()
```

## Value

A named list with the following elements:

- km_type - String. The type of k-means clustering. One of
  `c("minibatch", "standard")`. Defaults to `"minibatch"`.

- n_centroids - Integer or `NULL`. The number of centroids to use.
  `NULL` uses `sqrt(N_cells) * 4` centroids. Defaults to `NULL`.

- kmeans_iters - Integer. Number of maximum k-means iterations. Defaults
  to `100L`.

- batch_size - Integer. Max batch size, capped at `N_cells / 2`
  depending on the data set. Defaults to `4098L`.
