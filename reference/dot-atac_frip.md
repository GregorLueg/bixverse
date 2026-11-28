# Fraction of fragments in peaks from the ATAC library sizes

In the store's count unit: with insertion counts the share of insertions
in peaks (two per fragment), otherwise Signac's fragments in peaks over
fragments, which can exceed 1 when a fragment overlaps two peaks.

## Usage

``` r
.atac_frip(x, n_fragments)
```

## Arguments

- x:

  `SingleCellsAtac` or `SingleCellsMultiModal`.

- n_fragments:

  Numeric vector. Fragments per cell, every cell.

## Value

Numeric vector, or `NULL` for tiles, which cover everything.
