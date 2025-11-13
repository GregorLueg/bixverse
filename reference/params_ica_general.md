# Wrapper function for standard ICA parameters

Wrapper function for standard ICA parameters

## Usage

``` r
params_ica_general(maxit = 200L, alpha = 1, max_tol = 1e-04, verbose = FALSE)
```

## Arguments

- maxit:

  Integer. Maximum number of iterations for ICA.

- alpha:

  Float. The alpha parameter for the logcosh version of ICA. Should be
  between 1 to 2.

- max_tol:

  Numeric. Should be `0 < max_tol < 1`. Maximum tolerance of the
  algorithm.

- verbose:

  Boolean. Controls verbosity of the function.

## Value

A list with the parameters for usage in subsequent functions.
