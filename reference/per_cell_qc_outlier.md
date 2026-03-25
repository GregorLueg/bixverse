# Use MAD outlier detection on per-cell QC metrics

Use MAD outlier detection on per-cell QC metrics

## Usage

``` r
per_cell_qc_outlier(
  metric,
  threshold = 3,
  direction = c("twosided", "below", "above")
)
```

## Arguments

- metric:

  Numerical vector. The QC metric to check for.

- threshold:

  Numeric. How many MADs in either direction to consider for outlier
  detection.

- direction:

  String. One of `c("twosided", "below", "above")`. Which directionality
  to consider

## Value

A list with:

- outlier - Boolean vector indicating which cell is an outlier

- metrics - The applied thresholds.
