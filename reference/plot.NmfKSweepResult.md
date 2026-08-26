# Plot a consensus NMF k sweep

Stability and relative reconstruction error against `k`, the standard
cNMF diagnostic. Pick the `k` where stability is still high and the
error curve has not yet flattened out. Values of `k` whose consensus
step failed have no stability point.

## Usage

``` r
# S3 method for class 'NmfKSweepResult'
plot(x, ...)
```

## Arguments

- x:

  `NmfKSweepResult` object.

- ...:

  Additional params. Currently unused.

## Value

A `ggplot2` object with the two curves side by side.

## References

Kotliar et al., eLife, 2019
