# Calculate the critical value

This function calculates the critical value for a given set based on
random permutations and a given alpha value.

## Usage

``` r
rs_critval(values, iters, alpha, seed)
```

## Arguments

- values:

  Numeric vector. The full data set for which to calculate the critical
  value.

- iters:

  Integer. Number of random permutations to use.

- alpha:

  Float. The alpha value. For example, 0.001 would mean that the
  critical value is smaller than 0.1 percentile of the random
  permutations.

- seed:

  Integer. For reproducibility purposes

## Value

The critical value for the given parameters.
