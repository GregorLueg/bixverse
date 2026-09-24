# Calculate the column wise differential correlation between two sets of data.

**\[experimental\]** This function calculates the differential
correlation based on the Fisher method. For speed purposes, the function
will only calculate the differential correlation on the upper triangle
of the two correlation matrices.

## Usage

``` r
rs_differential_cor(x_a, x_b, spearman)
```

## Arguments

- x_a:

  Numeric matrix a, samples x features.

- x_b:

  Numeric matrix b, samples x features. Needs the same number of columns
  as `x_a`.

- spearman:

  Boolean. Shall the Spearman correlation be calculated instead of
  Pearson.

## Value

A list containing, one entry per upper-triangle feature pair:

- r_a - The correlation coefficients in the upper triangle of matrix a.

- r_b - The correlation coefficients in the upper triangle of matrix b.

- z_score - The z-scores of the difference in correlation coefficients.

- p_val - The z-scores transformed to two-sided p-values.
