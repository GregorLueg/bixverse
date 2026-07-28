# Run multiple NMF (HALS) restarts on a bulk expression matrix

**\[experimental\]** Runs `n_runs` HALS-NMF fits with random
initialisations seeded by `seed + i`. Expects samples x features.

## Usage

``` r
rs_nmf_multi_bulk(x, k, preprocessing, nmf_hals_params, n_runs, seed, verbose)
```

## Arguments

- x:

  Numerical matrix. Rows = samples, columns = features.

- k:

  Integer. Number of latent factors per run.

- preprocessing:

  String. One of `c("none", "sd", "sqrt_sd")`.

- nmf_hals_params:

  Named list. See
  [`params_nmf_hals()`](https://gregorlueg.github.io/bixverse/reference/params_nmf_hals.md).

- n_runs:

  Integer. Number of random restarts.

- seed:

  Integer. Base random seed. Run `i` uses `seed + i`.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with the following items

- w_all - Column-bound `W` matrices across all runs, shape
  `n_samples x (k * n_runs)`.

- h_per_run - List of `H` matrices, each `k x n_features`.

- losses - Numeric vector. Final reconstruction loss per run.

- converged - Logical vector. Convergence flag per run.

- best_idx - Integer. 1-indexed position of the run with the lowest
  final loss.
