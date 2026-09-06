# Prepare class for ICA

This is the generic function for doing the necessary preprocessing for
running independent component analysis.

## Usage

``` r
ica_processing(object, fast_svd = TRUE, random_seed = 123L, .verbose = TRUE)
```

## Arguments

- object:

  The class, see
  [`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md).
  Ideally, you should run
  [`preprocess_bulk_coexp()`](https://gregorlueg.github.io/bixverse/reference/preprocess_bulk_coexp.md)
  before applying this function.

- fast_svd:

  Boolean. Shall randomised SVD be used for the whitening. This is
  faster and usually causes little precision loss.

- random_seed:

  Integer. Seed for the randomised SVD. Only relevant, if fast_svd =
  `TRUE`.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

`BulkCoExp` with the needed data for ICA in the properties of the class.

## Examples

``` r
# whitening ahead of the ICA runs
mat <- t(synthetic_signal_matrix()$mat)
obj <- BulkCoExp(mat, data.table::data.table(sample_id = rownames(mat)))
obj <- preprocess_bulk_coexp(obj, hvg = 0.3, .verbose = FALSE)
obj <- ica_processing(obj, .verbose = FALSE)
dim(obj@processed_data$K)
#> [1]  89 300
```
