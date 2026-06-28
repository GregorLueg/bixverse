# K-mean parameter defaults.

Helper function to generate defaults for the k-mean clustering were more
control is needed.

## Usage

``` r
params_kmeans_defaults()
```

## Value

A list with the following parameters

- k_means_iter - Integer. The number of iterations to use for the
  clustering.

- k_means_init - String. The initialisation. Options are `"random"` and
  `"parallel"`. Defaults to `"parallel"`.

- gemm - Optional boolean. Controls which CPU implementation is used by
  the method. GEMM is faster with large dimensionality.

- hamerly - Optional boolean. Shall a faster exact method be used
  leveraging the triangle inequality. Faster on large data sets with
  large numbers of centroids.
