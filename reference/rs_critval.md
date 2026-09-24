# Calculate the critical value

**\[experimental\]** This function calculates the critical value for a
given set based on a bootstrap sample (with replacement) and a given
alpha value. The sample is sorted in decreasing order and the value at
the upper `alpha` tail is returned.

## Usage

``` r
rs_critval(values, iters, alpha, seed)
```

## Arguments

- values:

  Numeric vector. The full data set for which to calculate the critical
  value.

- iters:

  Integer. Size of the bootstrap sample.

- alpha:

  Float. The alpha value. For example, 0.001 would return the value
  exceeded by roughly 0.1 percent of the bootstrap sample.

- seed:

  Integer. For reproducibility purposes

## Value

The critical value for the given parameters.
