# Wrapper function for standard ICA parameters

Wrapper function for standard ICA parameters

## Usage

``` r
params_ica_general(maxit = 200L, alpha = 1, max_tol = 1e-04, verbose = FALSE)
```

## Arguments

- maxit:

  Integer. Maximum number of iterations for ICA. Defaults to `200L`.

- alpha:

  Numeric. The alpha parameter for the logcosh version of ICA. Should be
  between 1 to 2. Defaults to `1.0`.

- max_tol:

  Numeric. Should be `0 < max_tol < 1`. Maximum tolerance of the
  algorithm. Defaults to `1e-04`.

- verbose:

  Boolean. Controls verbosity of the function. Defaults to `FALSE`.

## Value

A named list with the following elements:

- maxit - Integer. Maximum number of iterations for ICA. Defaults to
  `200L`.

- alpha - Numeric. The alpha parameter for the logcosh version of ICA.
  Should be between 1 to 2. Defaults to `1.0`.

- max_tol - Numeric. Should be `0 < max_tol < 1`. Maximum tolerance of
  the algorithm. Defaults to `1e-04`.

- verbose - Boolean. Controls verbosity of the function. Defaults to
  `FALSE`.
