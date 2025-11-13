# Run ICA with cross-validation and random initialsiation

This function will split the data into `no_folds` and apply ICA with
`no_random_inits` over that fold.

## Usage

``` r
rs_ica_iters_cv(
  x,
  no_comp,
  no_folds,
  no_random_init,
  ica_type,
  random_seed,
  ica_params
)
```

## Arguments

- x:

  Numeric matrix. The processed data (no whitening function has been
  applied yet.)

- no_comp:

  Integer. Number of components to test for.

- no_folds:

  Integer. Number of folds to use for the cross-validation.

- no_random_init:

  Integer. Number of random initialisations.

- ica_type:

  String. Which type of ICA shall be run.

- random_seed:

  Integer. For reproducibility.

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
