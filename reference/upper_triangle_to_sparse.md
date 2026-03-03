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

  Integer. Did you shift the diagonal up.

- n:

  Integer. Number of columns and rows of the symmetric matrix.

- type:

  String. One of `c("csc", "csr")`. Which type of of sparse matrix to
  return.

## Value

The sparse matrix.
