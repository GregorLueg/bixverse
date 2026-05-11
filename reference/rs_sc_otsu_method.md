# Run Otsu's method

Maximises between-class variance of the observed score distribution to
find the optimal binary split. Robust to both bimodal and skewed
distributions.

## Usage

``` r
rs_sc_otsu_method(scores, bins)
```

## Arguments

- scores:

  Numeric vector. The vector for which to identify the threshold.

- bins:

  Integer. Number of bins to use for the histogram building.

## Value

The threshold based on Otsu's method
