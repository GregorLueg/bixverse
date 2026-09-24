# Run randomised SVD over a matrix

**\[experimental\]** Runs a randomised singular value decomposition over
a matrix. This implementation is faster than the full SVD on large data
sets, with slight loss in precision.

## Usage

``` r
rs_random_svd(x, scale, rank, seed, oversampling, n_power_iter)
```

## Arguments

- x:

  Numeric matrix. Rows = samples, columns = features.

- scale:

  Boolean. Shall the columns be variance normalised. (Mean centring will
  automatically occur.)

- rank:

  Integer. The rank to use.

- seed:

  Integer. Random seed for reproducibility.

- oversampling:

  Optional integer. Defaults to `10L` if `NULL`.

- n_power_iter:

  Optional integer. Number of power iterations (each with a QR
  decomposition). Defaults to `2L` if `NULL`.

## Value

A list with:

- scores - u matrix of the SVD multiplied by the singular values.

- v - v matrix of the SVD.

- s - Singular values of the SVD.

- scaled - Boolean. Was the matrix scaled.
