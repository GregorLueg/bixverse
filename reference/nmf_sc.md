# Run single-run NMF on single cell or meta cell data

Runs a single HALS NMF on a chosen subset of cells and genes. For
`SingleCells`, the counts are streamed from disk via the Rust binary
files; for `MetaCells`, the in-memory sparse counts are used.

## Usage

``` r
nmf_sc(
  object,
  k,
  cell_ids = NULL,
  gene_ids = NULL,
  preprocessing = "none",
  use_second_layer = TRUE,
  nmf_hals_params = params_nmf_hals(),
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

- seed:

  Integer. Random seed for initialisation.

- .verbose:

  Boolean or integer. Verbosity.

## Value

An `NmfResult` object.
