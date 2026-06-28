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

  Boolean. Shall the data be mean centered

- normalise_variance:

  Boolean. Shall the data have normalised variance

- randomised:

  Boolean. Shall fast, approximate randomised SVD be used.

- clr:

  Boolean. Shall the CLR-type `PFlogPF` be applied, see Booeshaghi, et
  al.

- size_factor:

  Numeric. The used size factor during I/O. It needs to be the same as
  during I/O to have correct results when using the `PFlogPF`
  transformation.

## Value

A list with the parameters
