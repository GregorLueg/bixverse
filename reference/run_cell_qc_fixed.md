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
