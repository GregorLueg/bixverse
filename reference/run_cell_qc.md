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
  mad = TRUE,
  r_standard = TRUE
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

- r_standard:

  Logical. If `TRUE`, the MAD is scaled by 1.4826 like
  [`stats::mad()`](https://rdrr.io/r/stats/mad.html). Set to `FALSE` to
  use the raw MAD (the previous behaviour).

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
#> CellQc: 100 cells, 2 outliers (2.0%)
#> Metrics:
#>   - lib_size: 2 outliers (mad = 2)
#>     MAD lower = 722.98
#>   - pct_mt: 0 outliers (mad = 0)
#>     MAD upper = 30.72
```
