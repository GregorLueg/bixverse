# Fixed-threshold cell QC

Thin wrapper around `run_cell_qc` with MAD disabled.

## Usage

``` r
run_cell_qc_fixed(metrics, cells_to_keep, hard_thresholds, groups = NULL)
```

## Arguments

- metrics:

  Named list of numeric vectors.

- cells_to_keep:

  Integer. 0-indexed cell positions.

- hard_thresholds:

  Required. See `run_cell_qc`.

- groups:

  Optional grouping vector.

## Examples

``` r
# a hard upper bound, no MAD anywhere
set.seed(42L)
qc <- run_cell_qc_fixed(
  metrics = list(pct_mt = runif(100, 0, 30)),
  cells_to_keep = 0:99,
  hard_thresholds = list(pct_mt = c(upper = 15))
)
sum(qc$combined)
#> [1] 55
```
