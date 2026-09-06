# Robust scaler.

Robust scaling, i.e., removes the median and scales data based on the
interquartile range (IQR). Useful if outliers are expected. NAs will be
ignored.

## Usage

``` r
robust_scale(x)
```

## Arguments

- x:

  Numeric vector.

## Value

x, robustly scaled.

## Examples

``` r
# median-centred, IQR-scaled vector
set.seed(123)
head(robust_scale(rnorm(10)))
#> [1] -0.5283040 -0.1652517  1.8010294  0.1652517  0.2298599  1.9728912
```
