# Generate a vector-based representation of the upper triangle of a matrix

**\[experimental\]** This function generates a vector from the upper
triangle of a given symmetric matrix. You have the option to remove the
diagonal with setting shift to 1.

## Usage

``` r
rs_dense_to_upper_triangle(x, shift)
```

## Arguments

- x:

  Numeric vector. The vector of correlation coefficients that you want
  to use to go back to a dense matrix.

- shift:

  Boolean. If you applied a shift, i.e. included the diagonal values. If
  `true`, assumes the diagonal values are `1`, otherwise derives them
  from the data.

## Value

The dense R matrix.
