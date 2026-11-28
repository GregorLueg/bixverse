# The ATAC feature table

The ATAC feature table

## Usage

``` r
.atac_var(x, cols = NULL)
```

## Arguments

- x:

  `SingleCellsAtac` or `SingleCellsMultiModal`.

- cols:

  Optional character vector of columns, in the `feature_idx` /
  `feature_id` naming.

## Value

A data.table keyed `feature_idx` / `feature_id` (1-based index).
