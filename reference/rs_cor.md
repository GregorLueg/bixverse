# Calculate the column wise correlations.

Calculates the correlation matrix of the columns. WARNING! Incorrect use
can cause kernel crashes. Wrapper around the Rust functions with type
checks are provided in the package.

## Usage

``` r
rs_cor(x, spearman)
```

## Arguments

- x:

  R matrix with doubles.

- spearman:

  Shall the Spearman correlation be calculated instead of Pearson.

## Value

The correlation matrix.
