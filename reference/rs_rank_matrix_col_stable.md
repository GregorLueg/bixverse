# Stable-gene rank matrix

**\[experimental\]** Unit-normalised column ranks computed against a set
of stable genes. Use with the `stable = TRUE` setting of the singscore
functions.

## Usage

``` r
rs_rank_matrix_col_stable(exp, stable_gene_indices)
```

## Arguments

- exp:

  Numerical matrix. The expression matrix (rows = genes, columns =
  samples).

- stable_gene_indices:

  Integer vector of stable genes. One-indexed.

## Value

Returns a matrix of normalised ranks with the same shape as `exp`.
