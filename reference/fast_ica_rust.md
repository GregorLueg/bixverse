# Fast ICA via Rust

This functions is a wrapper over the Rust implementation of fastICA and
the generation of the pre-processed data and pre-whitening matrix. It
has the same two options `c("logcosh", "exp")` to run ICA in parallel
modus. You can control the parameters of ICA via `ica_params`.

## Usage

``` r
fast_ica_rust(
  X,
  n_icas,
  ica_fun = c("logcosh", "exp"),
  ica_params = params_ica_general(),
  fast_svd = TRUE,
  seed = NULL
)
```

## Arguments

- X:

  Numeric matrix. The data on which you want to run fastICA.

- n_icas:

  Integer. Number of independent components to recover.

- ica_fun:

  String, element of `c("logcosh", "exp")`.

- ica_params:

  List. The ICA parameters, see
  [`params_ica_general()`](params_ica_general.md) wrapper function. This
  function generates a list containing:

  - maxit - Integer. Maximum number of iterations for ICA.

  - alpha - Float. The alpha parameter for the logcosh version of ICA.
    Should be between 1 to 2.

  - max_tol - Maximum tolerance of the algorithm.

  - verbose - Controls verbosity of the function.

- fast_svd:

  Boolean. Shall the randomised SVD be used. This is faster, but less
  precise.

- seed:

  Integer. Seed to ensure reproducible results.

## Value

A list containing:

- w The mixing matrix w.

- A ICA results matrix A.

- S ICA results matrix S.

- converged Boolean indicating if algorithm converged.
