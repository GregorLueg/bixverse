# Rust version of singscore for many gene sets

**\[experimental\]** Rust-based implementation of singscore over many
gene sets with optional paired down sets.

## Usage

``` r
rs_singscore_multi(
  ranks,
  up_list,
  down_list,
  center_score,
  known_direction,
  stable
)
```

## Arguments

- ranks:

  Numerical matrix. The ranked expression matrix.

- up_list:

  List. Up gene sets as zero-indexed indices. See
  [`rs_prepare_gsva_gs()`](https://gregorlueg.github.io/bixverse/reference/rs_prepare_gsva_gs.md).

- down_list:

  List or NULL. Paired down gene sets, same length and ordering as
  `up_list`.

- center_score:

  Boolean.

- known_direction:

  Boolean.

- stable:

  Boolean.

## Value

A named list with

- `scores` - Numerical matrix with the scores

- `dispersion` - Numerical matrix with the dispersions

Both matrices are of shape gene_sets × samples.
