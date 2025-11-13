# Run randomised SVD over a matrix

Runs a randomised singular value decomposition over a matrix. This
implementation is faster than the full SVD on large data sets, with
slight loss in precision.

## Usage

``` r
rs_random_svd(x, rank, seed, oversampling, n_power_iter)
```

## Arguments

- x:

  Numeric matrix. Rows = samples, columns = features.

- rank:

  Integer. The rank to use.

- seed:

  Integer. Random seed for reproducibility.

- oversampling:

  Integer. Defaults to `10L` if nothing is provided.

- n_power_iter:

  Integer. How often shall the QR decomposition be applied. Defaults to
  `2L` if nothing is provided.

## Value

A list with:

- u - u matrix of the SVD.

- v - v matrix of the SVD.

- s - Eigenvalues of the SVD.
