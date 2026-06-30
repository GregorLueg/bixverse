# Rust version of singscore for a single gene set

**\[experimental\]** Rust-based implementation of singscore for a single
up-regulated gene set with an optional paired down-regulated set.

## Usage

``` r
rs_singscore_single(
  ranks,
  up_set,
  down_set,
  center_score,
  known_direction,
  stable
)
```

## Arguments

- ranks:

  Numerical matrix. The ranked expression matrix (rows = genes, columns
  = samples). Produce with column-wise ranks (standard) or with
  [`rs_rank_matrix_col_stable()`](https://gregorlueg.github.io/bixverse/reference/rs_rank_matrix_col_stable.md)
  (stable).

- up_set:

  Integer vector. One-indexed gene indices of the up set.

- down_set:

  Integer vector or NULL. One-indexed gene indices of the optional down
  set.

- center_score:

  Boolean. Centre scores around 0. Disabled internally when
  `known_direction = FALSE`.

- known_direction:

  Boolean. Whether the up-set direction is known. Becomes irrelevant
  when `down_set` is also provided.

- stable:

  Boolean. If `TRUE`, use stable-gene score bounds.

## Value

A named list with `TotalScore`, `TotalDispersion`, and (when `down_set`
is provided) `UpScore`, `UpDispersion`, `DownScore`, `DownDispersion`.
