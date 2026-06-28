# Reconstruct a matrix from a flattened upper triangle vector

**\[experimental\]** This function takes a flattened vector of the upper
triangle from a symmetric matrix (think correlation matrix) and
reconstructs the full dense matrix for you.

## Usage

``` r
rs_upper_triangle_to_dense(data, shift, n)
```

## Arguments

- data:

  Numeric vector. The vector of for example correlation coefficients
  that you want to use to go back to a dense matrix.

- shift:

  Boolean. If you applied a shift, i.e. included the diagonal values. If
  `true`, assumes the diagonal values are `1`, otherwise derives them
  from the data.

- n:

  Integer. Original dimension (i.e., ncol/nrow) of the matrix to be
  reconstructed.

## Value

The dense R matrix.
