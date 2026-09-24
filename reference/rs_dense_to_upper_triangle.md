# Generate a vector-based representation of the upper triangle of a matrix

**\[experimental\]** This function generates a vector from the upper
triangle of a given symmetric matrix, iterating through the rows. You
have the option to remove the diagonal with setting `shift = TRUE`.

## Usage

``` r
rs_dense_to_upper_triangle(x, shift)
```

## Arguments

- x:

  Numeric matrix. The symmetric matrix to flatten.

- shift:

  Boolean. If `TRUE`, the diagonal is excluded, otherwise it is
  included.

## Value

Numeric vector with the upper triangle values.
