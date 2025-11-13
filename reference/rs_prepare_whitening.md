# Prepare the data for whitening

Prepares the data for subsequent usag in ICA. WARNING! Incorrect use can
cause kernel crashes. Wrapper around the Rust functions with type checks
are provided in the package.

## Usage

``` r
rs_prepare_whitening(x, fast_svd, seed, rank, oversampling, n_power_iter)
```

## Arguments

- x:

  The matrix to whiten. The whitening will happen over the columns.

- fast_svd:

  Boolean. Shall a randomised SVD be used. This is way faster on larger
  data sets.

- seed:

  Integer. Only relevant with fast_svd is set to `TRUE`.

- rank:

  Integer. How many ranks to use for the fast SVD approximation. If you
  supply `NULL`, it will default to `10L`. Only relevant with fast_svd
  is set to `TRUE`.

- oversampling:

  Integer. Oversampling parameter to make the approximation more
  precise. If you supply `NULL`, it will default to `10L`. Only relevant
  with fast_svd is set to `TRUE`.

- n_power_iter:

  Integer. How much shall the QR low rank approximation be powered. If
  you supply `NULL`, it will default to `2L`.

## Value

A list containing:

- x - The preprocessed matrix.

- k - The pre-whitening matrix k.
