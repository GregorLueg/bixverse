# Wrapper function for NMF (HALS) parameters

Wrapper function for NMF (HALS) parameters

## Usage

``` r
params_nmf_hals(
  max_iter = 250L,
  tol = 1e-04,
  eps = 1e-10,
  check_every = 10L,
  nmf_init = "nndsvd"
)
```

## Arguments

- max_iter:

  Integer. Maximum number of HALS iterations.

- tol:

  Numeric. Convergence tolerance on the relative change in
  reconstruction loss.

- eps:

  Numeric. Numerical floor for non-negativity / division safety.

- check_every:

  Integer. Convergence check interval in iterations.

- nmf_init:

  String. One of `c("nndsvd", "svd", "random")`. `"nndsvd"` and `"svd"`
  both map to deterministic NNDSVD initialisation; `"random"` uses
  random non-negative draws. For stabilised (multi-run) NMF this field
  is ignored and random init is always used.

## Value

A list with the HALS NMF parameters.
