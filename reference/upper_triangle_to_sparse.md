# Transform an upper triangle-stored matrix to a sparse one

Transform an upper triangle-stored matrix to a sparse one

## Usage

``` r
upper_triangle_to_sparse(upper_triangle_vals, shift, n, type = c("csc", "csr"))
```

## Arguments

- upper_triangle_vals:

  Numerical vector. The values of the upper triangle stored in a row
  major format.

- shift:

  Boolean. Did you exclude the diagonal.

- n:

  Integer. Number of columns and rows of the symmetric matrix.

- type:

  String. One of `c("csc", "csr")`. Which type of of sparse matrix to
  return.

## Value

The sparse matrix.

## Examples

``` r
# 3 x 3 symmetric matrix from its off-diagonal upper triangle
mat <- upper_triangle_to_sparse(c(0.5, 0.2, 0.8), shift = TRUE, n = 3L)
dim(mat)
#> [1] 3 3
```
