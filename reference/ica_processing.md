# Prepare class for ICA

This is the generic function for doing the necessary preprocessing for
running independent component analysis.

## Usage

``` r
ica_processing(object, fast_svd = TRUE, random_seed = 123L, .verbose = TRUE)
```

## Arguments

- object:

  The class, see [`bulk_coexp()`](bulk_coexp.md). Ideally, you should
  run [`preprocess_bulk_coexp()`](preprocess_bulk_coexp.md) before
  applying this function.

- fast_svd:

  Boolean. Shall randomised SVD be used for the whitening. This is
  faster and usually causes little precision loss.

- random_seed:

  Integer. Seed for the randomised SVD. Only relevant, if fast_svd =
  `TRUE`.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

`bulk_coexp` with the needed data for ICA in the properties of the
class.
