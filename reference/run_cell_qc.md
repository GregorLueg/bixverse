# Run MAD outlier detection on per-cell QC metrics

Run MAD outlier detection on per-cell QC metrics

## Usage

``` r
run_cell_qc(metrics, directions = NULL, threshold = 3)
```

## Arguments

- metrics:

  Named list of numeric vectors. Each element is a QC metric to check
  (e.g. `list(log10_lib_size = log10(lib_size), MT = mt_pct)`).

- directions:

  Named character vector mapping metric names to direction. One of
  `"twosided"`, `"below"`, `"above"`. Defaults to `"twosided"` for all
  metrics if `NULL`.

- threshold:

  Numeric. Number of MADs to use for outlier detection.

## Value

An object of class `CellQc` containing:

- metrics:

  The input metrics list.

- per_metric:

  Named list of per-metric results from
  [`per_cell_qc_outlier`](https://gregorlueg.github.io/bixverse/reference/per_cell_qc_outlier.md).

- outlier_mat:

  Logical matrix with one column per metric.

- combined:

  Logical vector. `TRUE` if a cell is an outlier in any metric.
