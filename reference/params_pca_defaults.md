# Helper function to generate default parameters for PCA

Helper function to generate default parameters for PCA

## Usage

``` r
params_pca_defaults()
```

## Value

A list with the following parameters for PCA.

- no_pcs - Integer. Number of PCs to consider. Defaults to `30L`.

- random_svd - Boolean. Shall randomised SVD be used. Defaults to
  `TRUE`.

- sparse - Boolean. Shall sparse solvers be used that do not do scaling.
  If set to yes, in the case of `random_svd = FALSE`, Lanczos iterations
  are used to solve the sparse SVD. With `random_svd = TRUE`, the sparse
  initial matrix is multiplied with the random matrix, yielding a much
  smaller dense matrix that does not increase the memory pressure
  massively.
