# Calculate the TOM from an expression matrix

Calculate the TOM from an expression matrix

## Usage

``` r
calculate_tom_from_exp(x, signed, version, cor_method)
```

## Arguments

- x:

  Numerical matrix. The expression matrix. Assumes that columns are the
  genes, and rows the samples.

- signed:

  Boolean. Do you want to calculate the signed version. If set to
  `FALSE`, the absolute correlation coefficients will be used.

- version:

  String. One of `c("v1", "v2")`. Defaults to `"v1"`

- cor_method:

  String. One of `c("pearson", spearman)`. Defaults to `"pearson"`.

## Value

The topological overlap matrix.

## Details

Calculates the topological overlap matrix from an expression matrix. It
will first calculate the specified correlation matrix and then generate
the TOM. The TOM is defined as:

**Unsigned, v1:**

\$\$TOM\_{ij} = \frac{a\_{ij} + \sum\_{k \neq i,j} a\_{ik}
a\_{kj}}{\min(k_i, k_j) + 1 - a\_{ij}}\$\$

**Signed, v1:**

\$\$TOM\_{ij} = \frac{a\_{ij} + \sum\_{k \neq i,j} a\_{ik}
a\_{kj}}{\min(k_i, k_j) + 1 - \left\|a\_{ij}\right\|}\$\$

**Unsigned, v2:**

\$\$TOM\_{ij} = 0.5 \left( a\_{ij} + \frac{\sum\_{k \neq i,j} a\_{ik}
a\_{kj}}{\min(k_i, k_j) + a\_{ij}} \right)\$\$

**Signed, v2:**

\$\$TOM\_{ij} = 0.5 \left( a\_{ij} + \frac{\sum\_{k \neq i,j} a\_{ik}
a\_{kj}}{\min(k_i, k_j) + \left\|a\_{ij}\right\|} \right)\$\$

where \\a\_{ij}\\ is the affinity between nodes \\i\\ and \\j\\, and
\\k_i = \sum_j a\_{ij}\\ is the connectivity of node \\i\\. For signed
networks, connectivity is calculated as \\k_i = \sum_j
\left\|a\_{ij}\right\|\\.

Version 2 uses a different normalization approach that scales the shared
neighbor contribution separately before combining it with the direct
connection strength.
