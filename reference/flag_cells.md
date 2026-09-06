# Add hard-threshold flags to a CellQc object

Unions new hard flags with existing ones per metric. Set `reset = TRUE`
to clear existing hard flags first.

## Usage

``` r
flag_cells(x, hard_thresholds, reset = FALSE)
```

## Arguments

- x:

  A `CellQc` object.

- hard_thresholds:

  Named list. See `run_cell_qc`.

- reset:

  Logical. Clear existing hard flags before applying.

## Value

Updated `CellQc`.

## Examples

``` r
# add a hard mitochondrial cut on top of the MAD flags
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
flag_cells(qc, list(pct_mt = c(upper = 15)))
#> CellQc: 100 cells, 25 outliers (25.0%)
#> Metrics:
#>   - lib_size: 5 outliers (mad = 5)
#>     MAD lower = 816.05
#>   - pct_mt: 21 outliers (mad = 0, hard = 21)
#>     MAD upper = 23.60
#>     Hard upper = 15.00
```
