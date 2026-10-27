# MAD outlier detection on per-group QC metrics

Aggregates each metric to per-group medians and applies
`per_cell_qc_outlier` to those medians, flagging whole groups (e.g.
donors, samples) whose median deviates by more than `threshold` MADs.
Directionality is respected per metric.

## Usage

``` r
per_group_qc_outlier(
  metrics,
  groups,
  directions,
  threshold = 3,
  r_standard = TRUE
)
```

## Arguments

- metrics:

  Named list of numeric vectors. The QC metrics.

- groups:

  Character vector. Group label per observation.

- directions:

  Named character vector mapping metric names to direction. One of
  `"twosided"`, `"below"`, `"above"`.

- threshold:

  Numeric. Number of MADs to use for outlier detection.

- r_standard:

  Boolean. Scale the MAD by 1.4826 like
  [`stats::mad()`](https://rdrr.io/r/stats/mad.html). See
  [`per_cell_qc_outlier()`](https://gregorlueg.github.io/bixverse/reference/per_cell_qc_outlier.md).

## Value

A `data.table` with one row per metric/group and columns `metric`,
`group`, `group_median`, `lower_threshold`, `upper_threshold`,
`is_outlier`.
