# Robust scaler.

Robust scaling, i.e., removes the median and scales data based on the
interquartile range (IQR). Useful if outliers are expected. NAs will be
ignored.

Robust scaling, i.e., removes the median and scales data based on the
interquartile range (IQR). Useful if outliers are expected. NAs will be
ignored.

## Usage

``` r
robust_scale(x)

robust_scale(x)
```

## Arguments

- x:

  Numeric vector.

## Value

x, robustly scaled.

x, robustly scaled.

## Examples

``` r
if (FALSE) { # \dontrun{

set.seed(123)
x <- rnorm(10)

x.scaled <- robust_scaling(x)
} # }
if (FALSE) { # \dontrun{

set.seed(123)
x <- rnorm(10)

x.scaled <- robust_scaling(x)
} # }
```
