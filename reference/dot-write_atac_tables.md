# Write the ATAC tables of a multi-modal object after a new store

Writes `var_atac`, the per-cell totals to obs and the count mode, resets
the ATAC state and cache, and recomputes a FRiP from an earlier fragment
QC against the new counts.

## Usage

``` r
.write_atac_tables(object, feature_ids, regions, lib_size, nnz, count_mode)
```

## Arguments

- object:

  `SingleCellsMultiModal` class.

- feature_ids:

  Character vector. One id per feature.

- regions:

  List with `chrom`, `start` and `end`.

- lib_size, nnz:

  Numeric vectors. Per-cell totals, every cell.

- count_mode:

  String. The count mode.

## Value

The object with the ATAC store attached.
