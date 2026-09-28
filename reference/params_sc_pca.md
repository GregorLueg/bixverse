# Wrapper for PCA specifically designed for single cells

Wrapper for PCA specifically designed for single cells

## Usage

``` r
params_sc_pca(
  mean_center = TRUE,
  normalise_variance = TRUE,
  randomised = NULL,
  clr = FALSE,
  size_factor = 10000,
  svd_solver = c("covariance", "randomised", "exact")
)
```

## Arguments

- mean_center:

  Boolean. Shall the data be mean centred Defaults to `TRUE`.

- normalise_variance:

  Boolean. Shall the data have normalised variance Defaults to `TRUE`.

- randomised:

  Boolean or `NULL`. Deprecated, use `svd_solver`. `TRUE` maps to
  `"randomised"`, `FALSE` to `"exact"`. Defaults to `NULL`.

- clr:

  Boolean. Shall the CLR-type `PFlogPF` be applied, see Booeshaghi, et
  al. Defaults to `FALSE`.

- size_factor:

  Numeric. The used size factor during I/O. It needs to be the same as
  during I/O to have correct results when using the `PFlogPF`
  transformation. Defaults to `10000.0`.

- svd_solver:

  String. Which solver to use. `"covariance"` builds the gene x gene
  cross-product and eigendecomposes it: exact, and the fastest option
  for a few thousand HVGs, but its cost grows with the square of the
  gene number in memory and the cube in time. `"randomised"` is a
  randomised SVD, approximate in the trailing components. `"exact"` is a
  full SVD on the dense path and Lanczos on the sparse one. One of
  `c("covariance", "randomised", "exact")`. Defaults to `"covariance"`.

## Value

A named list with the following elements:

- mean_center - Boolean. Shall the data be mean centred Defaults to
  `TRUE`.

- normalise_variance - Boolean. Shall the data have normalised variance
  Defaults to `TRUE`.

- randomised - Boolean or `NULL`. Deprecated, use `svd_solver`. `TRUE`
  maps to `"randomised"`, `FALSE` to `"exact"`. Defaults to `NULL`.

- clr - Boolean. Shall the CLR-type `PFlogPF` be applied, see
  Booeshaghi, et al. Defaults to `FALSE`.

- size_factor - Numeric. The used size factor during I/O. It needs to be
  the same as during I/O to have correct results when using the
  `PFlogPF` transformation. Defaults to `10000.0`.

- svd_solver - String. Which solver to use. `"covariance"` builds the
  gene x gene cross-product and eigendecomposes it: exact, and the
  fastest option for a few thousand HVGs, but its cost grows with the
  square of the gene number in memory and the cube in time.
  `"randomised"` is a randomised SVD, approximate in the trailing
  components. `"exact"` is a full SVD on the dense path and Lanczos on
  the sparse one. One of `c("covariance", "randomised", "exact")`.
  Defaults to `"covariance"`.
