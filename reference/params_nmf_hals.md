# Wrapper function for NMF (HALS) parameters

Wrapper function for NMF (HALS) parameters

## Usage

``` r
params_nmf_hals(
  max_iter = 250L,
  tol = 1e-04,
  eps = 1e-10,
  check_every = 10L,
  nmf_init = c("nndsvd", "svd", "random")
)
```

## Arguments

- max_iter:

  Integer. Maximum number of HALS iterations. Defaults to `250L`.

- tol:

  Numeric. Convergence tolerance on the relative change in
  reconstruction loss. Defaults to `1e-04`.

- eps:

  Numeric. Numerical floor for non-negativity / division safety.
  Defaults to `1e-10`.

- check_every:

  Integer. Convergence check interval in iterations. Defaults to `10L`.

- nmf_init:

  String. `"nndsvd"` and `"svd"` both map to deterministic NNDSVD
  initialisation; `"random"` uses random non-negative draws. For
  stabilised (multi-run) NMF this field is ignored and random init is
  always used. One of `c("nndsvd", "svd", "random")`. Defaults to
  `"nndsvd"`.

## Value

A named list with the following elements:

- max_iter - Integer. Maximum number of HALS iterations. Defaults to
  `250L`.

- tol - Numeric. Convergence tolerance on the relative change in
  reconstruction loss. Defaults to `1e-04`.

- eps - Numeric. Numerical floor for non-negativity / division safety.
  Defaults to `1e-10`.

- check_every - Integer. Convergence check interval in iterations.
  Defaults to `10L`.

- nmf_init - String. `"nndsvd"` and `"svd"` both map to deterministic
  NNDSVD initialisation; `"random"` uses random non-negative draws. For
  stabilised (multi-run) NMF this field is ignored and random init is
  always used. One of `c("nndsvd", "svd", "random")`. Defaults to
  `"nndsvd"`.
