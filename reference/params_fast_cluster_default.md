# Helper function to generate default parameters for the fast clustering for the doublet detection methods

Helper function to generate default parameters for the fast clustering
for the doublet detection methods

## Usage

``` r
params_fast_cluster_default()
```

## Value

A list with the following parameters for fast clustering:

- km_type - The type of k-means clustering. Defaults to `"minibatch"`

- n_centroids - The number of centroids to use. Default to `NULL` and
  the function will use `sqrt(N_cells) * 4` for the number of
  n_centroids.

- kmeans_iters - Number of maximum k-means iterations. Defaults to
  `100L`

- batch_size - Max batch size will be set to `4098L`, but pending data
  set set to `N_cells / 2`.
