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

## Examples

``` r
# MAD outlier detection over two metrics at once
set.seed(42L)
metrics <- list(
  lib_size = c(rnorm(99, 1000, 100), 50),
  pct_mt = runif(100, 0, 20)
)
run_cell_qc(
  metrics,
  cells_to_keep = 0:99,
  directions = c(lib_size = "below", pct_mt = "above")
)
#> CellQc: 100 cells, 5 outliers (5.0%)
#> Metrics:
#>   - lib_size: 5 outliers (mad = 5)
#>     MAD lower = 816.05
#>   - pct_mt: 0 outliers (mad = 0)
#>     MAD upper = 23.60
```
