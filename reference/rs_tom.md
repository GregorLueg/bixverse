# Calculates the TOM over an affinity matrix

Calculates the topological overlap measure for a given affinity matrix
x. Has the option to calculate the signed and unsigned version.

## Usage

``` r
rs_tom(x, tom_type, signed)
```

## Arguments

- x:

  Numerical matrix. Affinity matrix.

- tom_type:

  String. One of `c("v1", "v2")` - pending on choice, a different
  normalisation method will be used.

- signed:

  Boolean. Shall the signed TOM be calculated. If set to `FALSE`, values
  should be ≥ 0.

## Value

Returns the TOM matrix.
