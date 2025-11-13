# Calculate the column wise differential correlation between two sets of data.

This function calculates the differential correlation based on the
Fisher method. For speed purposes, the function will only calculate the
differential correlation on the upper triangle of the two correlation
matrices. WARNING! Incorrect use can cause kernel crashes. Wrapper
around the Rust functions with type checks are provided in the package.

## Usage

``` r
rs_differential_cor(x_a, x_b, spearman)
```

## Arguments

- x_a:

  R matrix a to be used for the differential correlation analysis.

- x_b:

  R matrix a to be used for the differential correlation analysis.

- spearman:

  Shall the Spearman correlation be calculated instead of Pearson.

## Value

A list containing:

- r_a - The correlation coefficients in the upper triangle of matrix a.

- r_b - The correlation coefficients in the upper triangle of matrix b.

- z_score - The z-scores of the difference in correlation coefficients.

- p_val - The z-scores transformed to p-values.
