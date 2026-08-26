# Plot the LDA topic count sweep

The three selection metrics and the combined score against `k`. Look for
the `k` where the combined score peaks, then sanity check that Cao Juan
has not already started climbing, which is the signal that topics are
duplicating.

## Usage

``` r
# S3 method for class 'LdaKSweepResult'
plot(x, ...)
```

## Arguments

- x:

  `LdaKSweepResult` object.

- ...:

  Ignored.

## Value

A `ggplot2` object with one panel per metric.
