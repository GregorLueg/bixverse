# Run multiple NMF (HALS) restarts on MetaCells

**\[experimental\]** Assumes that the sparse data is pre-filtered for
the cells/genes you wish to include. Indices in the sparse data need to
be 0-indexed.

## Usage

``` r
rs_nmf_multi_mc(
  sparse_data,
  k,
  preprocessing,
  use_second_layer,
  nmf_hals_params,
  n_runs,
  seed,
  verbose
)
```

## Arguments

- sparse_data:

  A named list with `data`, `indptr`, `indices`, `nrow`, `ncol` and
  `format`.

- k:

  Integer. Number of latent factors per run.

- preprocessing:

  String. One of `c("none", "sd", "sqrt_sd")`.

- use_second_layer:

  Boolean. If `TRUE`, runs NMF on normalised counts.

- nmf_hals_params:

  Named list. Contains the NMF parameters.

- n_runs:

  Integer. Number of random restarts.

- seed:

  Integer. Base random seed. Run `i` uses `seed + i`.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with `w_all`, `h_per_run`, `losses`, `converged`, `best_idx`
(1-indexed).
