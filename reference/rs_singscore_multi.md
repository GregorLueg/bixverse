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

  List. Up gene sets as 0-based indices. See
  [`rs_prepare_gsva_gs()`](https://gregorlueg.github.io/bixverse/reference/rs_prepare_gsva_gs.md).

- down_list:

  List or NULL. Paired down gene sets as 0-based indices, same length
  and ordering as `up_list`.

- center_score:

  Boolean. Centre scores around 0. Disabled internally when
  `known_direction = FALSE`.

- known_direction:

  Boolean. Whether the up-set direction is known.

- stable:

  Boolean. If `TRUE`, use stable-gene score bounds.

## Value

A named list with

- `scores` - Numerical matrix with the scores

- `dispersions` - Numerical matrix with the dispersions

Both matrices are of shape gene sets x samples.
