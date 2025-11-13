# Run the Rust implementation of fast ICA.

This function serves as a wrapper over the fast ICA implementations in
Rust. It assumes a whitened matrix and also an intialised w_init.
WARNING! Incorrect use can cause kernel crashes. Wrapper around the Rust
functions with type checks are provided in the package.

## Usage

``` r
rs_fast_ica(whiten, w_init, ica_type, ica_params)
```

## Arguments

- whiten:

  Numerical matrix. The whitened matrix.

- w_init:

  Numerical matrix. The initial unmixing matrix. ncols need to be equal
  to nrows of whiten.

- ica_type:

  String. One of 'logcosh' or 'exp'.

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

A list with the following items:

- mixing - The mixing matrix for subsequent usage.

- converged - Boolean if the algorithm converged.
