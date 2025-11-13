# Generate a vector-based representation of the upper triangle of a matrix

This function generates a vector from the upper triangle of a given
symmetric matrix. You have the option to remove the diagonal with
setting shift to 1.

## Usage

``` r
rs_dense_to_upper_triangle(x, shift)
```

## Arguments

- x:

  Numeric vector. The vector of correlation coefficients that you want
  to use to go back to a dense matrix.

- shift:

  Integer. If you want to apply a shift, i.e. included the diagonal
  values = 0; or excluded the diagonal values = 1.

## Value

The dense R matrix.
