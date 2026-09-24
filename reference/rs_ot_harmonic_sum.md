# Calculate the OT harmonic sum

**\[experimental\]** Sorts `x` in decreasing order and sums
`x[i] / i^2`, normalised by the same sum over a vector of ones of equal
length.

## Usage

``` r
rs_ot_harmonic_sum(x)
```

## Arguments

- x:

  The numeric vector (should be between 0 and 1) for which to calculate
  the harmonic sum

## Value

The normalised harmonic sum according to the OT calculation.
