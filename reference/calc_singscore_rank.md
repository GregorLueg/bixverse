# Rank an expression matrix for singscore

Column-ranks an expression matrix. If `stable_genes` is provided, uses
the stable-genes ranking method (unit-normalised ranks against a small
set of stable genes); otherwise standard column ranks. Sets
`attr(., "stable")` so downstream singscore functions know which bounds
formula to use.

## Usage

``` r
calc_singscore_rank(exp, stable_genes = NULL)
```

## Arguments

- exp:

  Numerical matrix. Rows = genes, columns = samples.

- stable_genes:

  Character vector or NULL. Gene names of stable genes. Defaults to
  `NULL` (standard ranking).

## Value

A rank matrix with the same shape as `exp` and `attr(., "stable")` set
to `TRUE` or `FALSE`.

## References

1.) Foroutan et al., BMC Bioinformatics, 2018.; 2.) Bhuva, et al.,
Nucleic Acids Res., 2020
