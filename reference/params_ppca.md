# Wrapper function for probabilistic PCA parameters

Parameters for PCA with missing values via probabilistic PCA, see
[`run_ppca()`](https://gregorlueg.github.io/bixverse/reference/run_ppca.md).
Defaults follow pcaMethods.

## Usage

``` r
params_ppca(
  n_pcs = 2L,
  max_iter = 1000L,
  tol = 1e-05,
  seed = 42L,
  centre = TRUE,
  scale = FALSE
)
```

## Arguments

- n_pcs:

  Integer. Number of principal components. Defaults to `2L`.

- max_iter:

  Integer. Maximum number of EM iterations. Defaults to `1000L`.

- tol:

  Numeric. Relative change in the objective below which EM stops.
  Defaults to `1e-05`.

- seed:

  Integer. Seed for the random initial loadings. Defaults to `42L`.

- centre:

  Boolean. Shall the observed column means be subtracted first. Defaults
  to `TRUE`.

- scale:

  Boolean. Shall the columns be divided by their observed standard
  deviation first (pcaMethods' `"uv"`). Defaults to `FALSE`.

## Value

A named list with the following elements:

- n_pcs - Integer. Number of principal components. Defaults to `2L`.

- max_iter - Integer. Maximum number of EM iterations. Defaults to
  `1000L`.

- tol - Numeric. Relative change in the objective below which EM stops.
  Defaults to `1e-05`.

- seed - Integer. Seed for the random initial loadings. Defaults to
  `42L`.

- centre - Boolean. Shall the observed column means be subtracted first.
  Defaults to `TRUE`.

- scale - Boolean. Shall the columns be divided by their observed
  standard deviation first (pcaMethods' `"uv"`). Defaults to `FALSE`.

## References

Tipping and Bishop, J R Stat Soc B, 1999; Stacklies, et al.,
Bioinformatics, 2007
