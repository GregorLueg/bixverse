# Find a threshold via the Otsu method

Find a threshold via the Otsu method

## Usage

``` r
find_threshold_otsu(x, bins = 100L)
```

## Arguments

- x:

  Numerical vector. The vector for which to find the threshold via the
  Otsu method

- bins:

  Integer. Number of bins to use.

## Value

The threshold.

## Examples

``` r
# threshold separating two Gaussian modes
set.seed(42)
find_threshold_otsu(c(rnorm(100), rnorm(100, mean = 5)))
#> [1] 2.300925
```
