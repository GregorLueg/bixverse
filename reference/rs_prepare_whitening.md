# Prepare the data for whitening

**\[experimental\]** Prepares the data for subsequent usage in ICA.
Incorrect use can cause kernel crashes. Wrapper around the Rust
functions with type checks are provided in the package.

## Usage

``` r
rs_prepare_whitening(x, fast_svd, seed, rank, oversampling, n_power_iter)
```

## Arguments

- x:

  Numeric matrix to whiten. The columns are centred and the whitening
  happens over the columns.

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

  Integer. Number of power iterations for the randomised SVD. If you
  supply `NULL`, it will default to `2L`. Only relevant with fast_svd is
  set to `TRUE`.

## Value

A list containing:

- x - The column-centred input, transposed.

- k - The whitening matrix K. With `fast_svd = TRUE` it can carry more
  than `rank` rows; the caller trims it.
