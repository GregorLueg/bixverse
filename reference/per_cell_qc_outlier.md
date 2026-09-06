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

## Examples

``` r
# one badly undersequenced cell, flagged on the lower tail only
set.seed(42L)
lib_size <- c(rnorm(99, 1000, 100), 50)
res <- per_cell_qc_outlier(lib_size, direction = "below")
sum(res$outlier)
#> [1] 5
res$metrics
#>          median upper_threshold lower_threshold 
#>       1008.9100              NA        816.0509 
```
