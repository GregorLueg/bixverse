# K-mean parameter defaults.

Helper function to generate defaults for the k-mean clustering were more
control is needed.

## Usage

``` r
params_kmeans_defaults()
```

## Value

A named list with the following elements:

- k_means_algorithm - String. The k-means algorithm. `"emeans"` (Bessa
  et al., 2026) runs mini-batch EMA updates and ignores `k_means_iter`,
  `k_means_init`, `gemm` and `hamerly`. `"lloyd"` is Lloyd's algorithm.
  One of `c("emeans", "lloyd")`. Defaults to `"emeans"`.

- e_means_epochs - Numeric or `NULL`. e-means only. Passes over the
  data. `NULL` takes the Rust default of 5. Defaults to `NULL`.

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
