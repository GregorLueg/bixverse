# Run NMF (HALS) on a bulk expression matrix (single run)

**\[experimental\]** Runs a single HALS-NMF fit on a dense matrix.
Expects samples x features. The resulting decomposition `V ~ W H` places
`W` (samples x k) as sample-side factors and `H` (k x features) as
feature loadings.

## Usage

``` r
rs_nmf_single_bulk(x, k, preprocessing, nmf_hals_params, seed, verbose)
```

## Arguments

- x:

  Numerical matrix. Rows = samples, columns = features.

- k:

  Integer. Number of latent factors.

- preprocessing:

  String. One of `c("none", "sd", "sqrt_sd")`.

- nmf_hals_params:

  Named list. See
  [`params_nmf_hals()`](https://gregorlueg.github.io/bixverse/reference/params_nmf_hals.md).

- seed:

  Integer. Random seed for the NMF initialisation.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with the following items

- w - The `W` matrix of shape `n_samples x k`.

- h - The `H` matrix of shape `k x n_features`.

- final_loss - Final reconstruction loss.

- n_iter - Number of iterations the algorithm ran for.

- converged - Did the NMF algorithm converge.
