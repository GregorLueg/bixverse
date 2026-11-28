# Write obs, var and the count mode of a fresh ATAC-only store

Write obs, var and the count mode of a fresh ATAC-only store

## Usage

``` r
.finish_atac_root(
  object,
  barcodes,
  feature_ids,
  regions,
  lib_size,
  nnz,
  count_mode,
  .verbose
)
```

## Arguments

- object:

  `SingleCellsAtac` class.

- barcodes:

  Character vector. Cells in store order.

- feature_ids:

  Character vector. One id per feature.

- regions:

  List with `chrom`, `start` and `end`.

- lib_size, nnz:

  Numeric vectors. Per-cell totals.

- count_mode:

  String. The count mode.

- .verbose:

  Boolean. Controls verbosity.

## Value

The loaded `SingleCellsAtac`.
