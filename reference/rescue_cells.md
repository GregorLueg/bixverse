# Rescue MAD-flagged cells that fall within safe bounds

For each metric in `rescue_thresholds`, cells whose value is within
`[lower, upper]` are un-flagged from MAD. Hard-threshold flags are not
cleared; a warning is issued when overlaps occur. Calling `rescue_cells`
overwrites any prior rescue set for the metrics named.

## Usage

``` r
rescue_cells(x, rescue_thresholds)
```

## Arguments

- x:

  A `CellQc` object.

- rescue_thresholds:

  Named list of numeric vectors with `lower` and/or `upper`.

## Value

Updated `CellQc`.

## Examples

``` r
# un-flag MAD outliers that still sit above a sane library size
set.seed(42L)
metrics <- list(
  lib_size = c(rnorm(99, 1000, 100), 50),
  pct_mt = runif(100, 0, 20)
)
qc <- run_cell_qc(
  metrics,
  cells_to_keep = 0:99,
  directions = c(lib_size = "below", pct_mt = "above")
)
rescue_cells(qc, list(lib_size = c(lower = 500)))
#> CellQc: 100 cells, 1 outliers (1.0%)
#> Metrics:
#>   - lib_size: 1 outliers (mad = 5, rescued = 4)
#>     MAD lower = 816.05
#>   - pct_mt: 0 outliers (mad = 0)
#>     MAD upper = 23.60
```
