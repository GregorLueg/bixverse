# Reconstruct a matrix from a flattened upper triangle vector

This function takes a flattened vector of the upper triangle from a
symmetric matrix (think correlation matrix) and reconstructs the full
dense matrix for you.

## Usage

``` r
rs_upper_triangle_to_dense(cor_vector, shift, n)
```

## Arguments

- cor_vector:

  Numeric vector. The vector of correlation coefficients that you want
  to use to go back to a dense matrix.

- shift:

  Integer. If you applied a shift, i.e. included the diagonal values =
  0; or excluded the diagonal values = 1.

- n:

  Integer. Original dimension (i.e., ncol/nrow) of the matrix to be
  reconstructed.

## Value

The dense R matrix.
