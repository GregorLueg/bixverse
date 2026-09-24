# K-mean parameter defaults.

Helper function to generate defaults for the k-mean clustering were more
control is needed.

## Usage

``` r
params_kmeans_defaults()
```

## Value

A named list with the following elements:

- k_means_iter - Integer. The number of iterations to use for the
  clustering. Defaults to `30L`.

- k_means_init - String. The initialisation. One of
  `c("parallel", "random")`. Defaults to `"parallel"`.

- gemm - Boolean or `NULL`. Controls which CPU implementation is used by
  the method. GEMM is faster with large dimensionality. Defaults to
  `FALSE`.

- hamerly - Boolean or `NULL`. Shall a faster exact method be used
  leveraging the triangle inequality. Faster on large data sets with
  large numbers of centroids. Defaults to `TRUE`.
