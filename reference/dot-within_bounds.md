# Check whether values fall within named bounds

Check whether values fall within named bounds

## Usage

``` r
.within_bounds(x, bounds)
```

## Arguments

- x:

  Numeric vector.

- bounds:

  Named numeric vector with `lower` and/or `upper`. At least one must be
  present.

## Value

Logical vector of length `length(x)`. `TRUE` where the value is within
the supplied bounds (inclusive on both ends).
