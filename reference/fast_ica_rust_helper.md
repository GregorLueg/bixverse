# Fast ICA via Rust from processed data

This functions is a wrapper over the Rust implementation of fastICA
expected already the pre-processed matrix and the pre-whitening matrix.
It has the same two options `c("logcosh", "exp")` to run ICA in parallel
modus. You can control the parameters of ICA via `ica_params`.

## Usage

``` r
fast_ica_rust_helper(
  X,
  K,
  n_icas,
  ica_fun = c("logcosh", "exp"),
  seed = NULL,
  ica_params = params_ica_general()
)
```

## Arguments

- X:

  Numeric matrix. Processed data. Output of
  [`rs_prepare_whitening()`](rs_prepare_whitening.md).

- K:

  The K matrix. Pre-whitening matrix. Output of
  [`rs_prepare_whitening()`](rs_prepare_whitening.md).

- n_icas:

  Integer. Number of independent components to recover.

- ica_fun:

  String, element of `c("logcosh", "exp")`.

- seed:

  Integer. Seed to ensure reproducible results.

- ica_params:

  List. The ICA parameters, see
  [`params_ica_general()`](params_ica_general.md) wrapper function. This
  function generates a list containing:

  - maxit - Integer. Maximum number of iterations for ICA.

  - alpha - Float. The alpha parameter for the logcosh version of ICA.
    Should be between 1 to 2.

  - max_tol - Maximum tolerance of the algorithm.

  - verbose - Controls verbosity of the function.

## Value

A list containing:

- w The mixing matrix w.

- A ICA results matrix A.

- S ICA results matrix S.

- converged Boolean indicating if algorithm converged.
