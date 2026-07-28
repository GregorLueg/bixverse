# Run outlier detection on per-cell QC metrics

Helper function to run initial quality control on cells.

## Usage

``` r
run_cell_qc(
  metrics,
  cells_to_keep,
  directions = NULL,
  threshold = 3,
  groups = NULL,
  hard_thresholds = NULL,
  mad = TRUE
)
```

## Arguments

- metrics:

  Named list of numeric vectors.

- cells_to_keep:

  Integer. 0-indexed cell positions.

- directions:

  Named character vector, one of `"twosided"`, `"below"`, `"above"`.
  Defaults to `"twosided"`.

- threshold:

  Numeric. MADs for MAD outlier detection.

- groups:

  Optional grouping vector.

- hard_thresholds:

  Optional named list of numeric vectors with `lower` and/or `upper`
  bounds, e.g. `list(MT = c(upper = 15))`. Applied independent of
  groups.

- mad:

  Logical. If `FALSE`, skip MAD entirely; `hard_thresholds` must then be
  supplied.

## Value

A `CellQc` object.
