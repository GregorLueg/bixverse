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
