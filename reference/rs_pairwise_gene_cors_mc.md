# Calculate the pairwise gene-correlation for meta cells

**\[experimental\]**

## Usage

``` r
rs_pairwise_gene_cors_mc(sparse_data, gene_indices_1, gene_indices_2, spearman)
```

## Arguments

- sparse_data:

  A named list that needs to have `data`, `indptr`, `indices`, `nrow`,
  `ncol` and `format`.

- gene_indices_1:

  Integer. The gene indices for the first set of genes. Must be
  0-indexed!

- gene_indices_2:

  Integer. The gene indices for the first set of genes. Must be
  0-indexed!

- spearman:

  Boolean. Shall the spearman correlation be calculated.

## Value

The vector of correlations between the pairs of gene_indices_1 and
gene_indices_2
