# Run ICA over a given no_comp with random initilisations of w_init

This function implements a stabilised ICA like algorithm in Rust.
Briefly, it generates random w_init matrices (total number being
no_random_init) and runs ICA given the x_processed and k data over
these. The function returns combined S from the different runs and a
boolean vector indicating if this specific run converged.

## Usage

``` r
rs_ica_iters(x1, k, no_comp, no_random_init, ica_type, random_seed, ica_params)
```

## Arguments

- x1:

  Numerical matrix. The processed matrix (but not yet whitened!)

- k:

  Numerical matrix. The whitening matrix.

- no_comp:

  Integer. Number of independent components to return.

- no_random_init:

  Integer. Number of random initialisations to test.

- ica_type:

  String. One of 'logcosh' or 'exp'.

- random_seed:

  Integer. Seed for randomisations.

- ica_params:

  A list containing:

  - maxit - Integer. Maximum number of iterations for ICA.

  - alpha - Float. The alpha parameter for the logcosh version of ICA.
    Should be between 1 to 2.

  - max_tol - Maximum tolerance of the algorithm

  - verbose - Verbosity of the function, i.e., shall individual iters be
    shown.

  If the list is empty or the expected elements are not found, default
  values are used.

## Value

A list containing:

- s_combined - The combined matrices for S. Dimensions are nrows =
  features; and ncols = ncomp \* no_random_init.

- converged - Boolean vector indicating if the respective run reached
  convergence. Length = no_random_init
