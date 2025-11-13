# Calculate the column wise correlations.

Calculates the correlation matrix of the columns. This function will
return the upper triangle. WARNING! Incorrect use can cause kernel
crashes. Wrapper around the Rust functions with type checks are provided
in the package.

## Usage

``` r
rs_cor_upper_triangle(x, spearman, shift)
```

## Arguments

- x:

  R matrix with doubles.

- spearman:

  Shall the Spearman correlation be calculated instead of Pearson.

- shift:

  Shall a shift be applied to the matrix. 0 = the diagonal will be
  included. 1 = the diagonal will not be included.

## Value

The upper triangle of the correlation matrix iterating through the rows,
shifted by one (the diagonal will not be returned).
