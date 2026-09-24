# Wrapper for PCA specifically designed for single cells

Wrapper for PCA specifically designed for single cells

## Usage

``` r
params_sc_pca(
  mean_center = TRUE,
  normalise_variance = TRUE,
  randomised = TRUE,
  clr = FALSE,
  size_factor = 10000
)
```

## Arguments

- mean_center:

  Boolean. Shall the data be mean centred Defaults to `TRUE`.

- normalise_variance:

  Boolean. Shall the data have normalised variance Defaults to `TRUE`.

- randomised:

  Boolean. Shall fast, approximate randomised SVD be used. Defaults to
  `TRUE`.

- clr:

  Boolean. Shall the CLR-type `PFlogPF` be applied, see Booeshaghi, et
  al. Defaults to `FALSE`.

- size_factor:

  Numeric. The used size factor during I/O. It needs to be the same as
  during I/O to have correct results when using the `PFlogPF`
  transformation. Defaults to `10000.0`.

## Value

A named list with the following elements:

- mean_center - Boolean. Shall the data be mean centred Defaults to
  `TRUE`.

- normalise_variance - Boolean. Shall the data have normalised variance
  Defaults to `TRUE`.

- randomised - Boolean. Shall fast, approximate randomised SVD be used.
  Defaults to `TRUE`.

- clr - Boolean. Shall the CLR-type `PFlogPF` be applied, see
  Booeshaghi, et al. Defaults to `FALSE`.

- size_factor - Numeric. The used size factor during I/O. It needs to be
  the same as during I/O to have correct results when using the
  `PFlogPF` transformation. Defaults to `10000.0`.
