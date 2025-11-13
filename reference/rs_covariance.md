# Calculate the column-wise co-variance.

Calculates the co-variance of the columns. WARNING! Incorrect use can
cause kernel crashes. Wrapper around the Rust functions with type checks
are provided in the package.

## Usage

``` r
rs_covariance(x)
```

## Arguments

- x:

  R matrix with doubles.

## Value

The co-variance matrix.
