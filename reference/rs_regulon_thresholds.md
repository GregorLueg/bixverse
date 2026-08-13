# Derive the on/off threshold per regulon in Rust

**\[experimental\]** Fits a two-component Gaussian mixture per regulon
and compares it against a single Gaussian by BIC. If the mixture wins,
the threshold is the kernel density minimum between the two component
means, otherwise it falls back to `mean + 2 * sd`. This follows pySCENIC
rather than AUCell.

## Usage

``` r
rs_regulon_thresholds(auc_matrix, binarise_params)
```

## Arguments

- auc_matrix:

  Numeric matrix of cells x regulons, as returned by
  [`rs_aucell()`](https://gregorlueg.github.io/bixverse/reference/rs_aucell.md).

- binarise_params:

  List. The binarisation parameters, see
  [`params_scenic_binarise()`](https://gregorlueg.github.io/bixverse/reference/params_scenic_binarise.md).

## Value

A list with `thresholds` (one per regulon) and `bimodal` (whether the
mixture won the BIC comparison).
