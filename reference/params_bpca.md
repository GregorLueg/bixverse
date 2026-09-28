# Wrapper function for Bayesian PCA parameters

Parameters for PCA with missing values via Bayesian PCA, see
[`run_bpca()`](https://gregorlueg.github.io/bixverse/reference/run_bpca.md).
Defaults follow pcaMethods.

## Usage

``` r
params_bpca(
  n_pcs = 2L,
  max_iter = 100L,
  tol = 1e-04,
  centre = TRUE,
  scale = FALSE
)
```

## Arguments

- n_pcs:

  Integer. Number of principal components. Defaults to `2L`.

- max_iter:

  Integer. Maximum number of variational steps. Defaults to `100L`.

- tol:

  Numeric. Change in `log10(tau)` over ten steps below which the fit
  stops. Defaults to `1e-04`.

- centre:

  Boolean. Shall the observed column means be subtracted first. Defaults
  to `TRUE`.

- scale:

  Boolean. Shall the columns be divided by their observed standard
  deviation first (pcaMethods' `"uv"`). Defaults to `FALSE`.

## Value

A named list with the following elements:

- n_pcs - Integer. Number of principal components. Defaults to `2L`.

- max_iter - Integer. Maximum number of variational steps. Defaults to
  `100L`.

- tol - Numeric. Change in `log10(tau)` over ten steps below which the
  fit stops. Defaults to `1e-04`.

- centre - Boolean. Shall the observed column means be subtracted first.
  Defaults to `TRUE`.

- scale - Boolean. Shall the columns be divided by their observed
  standard deviation first (pcaMethods' `"uv"`). Defaults to `FALSE`.

## References

Oba, et al., Bioinformatics, 2003; Stacklies, et al., Bioinformatics,
2007
