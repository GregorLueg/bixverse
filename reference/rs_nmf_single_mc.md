# Run NMF (HALS) on MetaCells

**\[experimental\]** Assumes that the sparse data is pre-filtered for
the cells/genes you wish to include. Indices in the sparse data need to
be 0-indexed.

## Usage

``` r
rs_nmf_single_mc(
  sparse_data,
  k,
  preprocessing,
  use_second_layer,
  nmf_hals_params,
  seed,
  verbose
)
```

## Arguments

- sparse_data:

  A named list with `data`, `indptr`, `indices`, `nrow`, `ncol` and
  `format`.

- k:

  Integer. Number of latent factors to return.

- preprocessing:

  String. One of `c("none", "sd", "sqrt_sd")`.

- use_second_layer:

  Boolean. If `TRUE`, runs NMF on normalised counts.

- nmf_hals_params:

  Named list. Contains the NMF parameters.

- seed:

  Integer. Random seed for initialisation.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with `w`, `h`, `final_loss`, `n_iter`, `converged`.
