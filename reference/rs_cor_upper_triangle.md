# Calculate the column wise correlations and returns the upper triangle

**\[experimental\]** Calculates the correlation matrix of the columns,
but returns the upper triangle only as a flat vector.

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

  Boolean. If `TRUE`, the diagonal is excluded, otherwise it is
  included.

## Value

The upper triangle of the correlation matrix iterating through the rows,
with or without the diagonal depending on `shift`.
