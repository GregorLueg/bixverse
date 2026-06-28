# Run MAD outlier detection on per-cell QC metrics

Run MAD outlier detection on per-cell QC metrics

## Usage

``` r
run_cell_qc(
  metrics,
  cells_to_keep,
  directions = NULL,
  threshold = 3,
  groups = NULL
)
```

## Arguments

- metrics:

  Named list of numeric vectors. Each element is a QC metric to check
  (e.g. `list(log10_lib_size = log10(lib_size), MT = mt_pct)`).

- cells_to_keep:

  Integer. Which cells were included in the analysis. 0-indices for
  Rust.

- directions:

  Named character vector mapping metric names to direction. One of
  `"twosided"`, `"below"`, `"above"`. Defaults to `"twosided"` for all
  metrics if `NULL`.

- threshold:

  Numeric. Number of MADs to use for outlier detection.

- groups:

  Optional grouping variable. An atomic vector of length equal to the
  metrics. Per-group outlier detection only runs when more than one
  group is present.

## Value

An object of class `CellQc` containing:

- cell_idx:

  Integer vector of 1-indexed cell positions.

- metrics:

  The input metrics list.

- groups:

  Character vector of group labels.

- per_metric:

  Named list of per-metric results from
  [`per_cell_qc_outlier`](https://gregorlueg.github.io/bixverse/reference/per_cell_qc_outlier.md).

- outlier_mat:

  Logical matrix with one column per metric.

- combined:

  Logical vector. `TRUE` if a cell is an outlier in any metric.

- per_group_stats:

  A `data.table` of per-group outlier statistics (see
  [`per_group_qc_outlier`](https://gregorlueg.github.io/bixverse/reference/per_group_qc_outlier.md)),
  or `NULL` for a single group.
