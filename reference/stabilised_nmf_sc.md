# Run stabilised (multi-run) NMF on single cell or meta cell data

Runs `n_runs` HALS NMF with random initialisations seeded by `seed + i`.
The `nmf_init` field in `nmf_hals_params` is ignored; random init is
always used.

## Usage

``` r
stabilised_nmf_sc(
  object,
  k,
  cell_ids = NULL,
  gene_ids = NULL,
  preprocessing = "none",
  use_second_layer = TRUE,
  nmf_hals_params = params_nmf_hals(),
  n_runs = 30L,
  seed = 42L,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells` or `MetaCells` class.

- k:

  Integer. Number of latent factors to return.

- cell_ids:

  Optional character. Cell ids (or meta cell ids) to restrict the NMF
  to. If `NULL`, uses
  [`get_cells_to_keep()`](https://gregorlueg.github.io/bixverse/reference/get_cells_to_keep.md)
  for `SingleCells` and all meta cells for `MetaCells`.

- gene_ids:

  Optional character. Gene ids to restrict the NMF to. If `NULL`, uses
  [`get_hvg()`](https://gregorlueg.github.io/bixverse/reference/get_hvg.md)
  on the object.

- preprocessing:

  String. One of `c("none", "sd", "sqrt_sd")`.

- use_second_layer:

  Boolean. If `TRUE`, runs NMF on the normalised counts (recommended);
  if `FALSE`, on the raw counts.

- nmf_hals_params:

  List, see
  [`params_nmf_hals()`](https://gregorlueg.github.io/bixverse/reference/params_nmf_hals.md).

- n_runs:

  Integer. Number of random restarts.

- seed:

  Integer. Random seed for initialisation.

- .verbose:

  Boolean or integer. Verbosity.

## Value

A `StabilisedNmfResult` object.
