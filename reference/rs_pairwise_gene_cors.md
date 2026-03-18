# Calculates pairwise gene correlations in single cell

Calculates pairwise gene correlations in single cell

## Usage

``` r
rs_pairwise_gene_cors(
  f_path,
  gene_indices_1,
  gene_indices_2,
  cells_to_keep,
  spearman
)
```

## Arguments

- f_path:

  Path to the `counts_genes.bin` file.

- gene_indices_1:

  Integer vector. The first set of gene indices to correlate against
  `gene_indices_2` (0-indexed!)

- gene_indices_2:

  Integer vector. The second set of gene indices. Needs to be of same
  length as `gene_indices_1`. (0-indexed!)

- cells_to_keep:

  Integer. The indices of the cells to include in this analysis.
  (0-indexed!)

- spearman:

  Boolean. Shall spearman correlation be used.

## Value

The vector of correlations between the pairs of gene_indices_1 and
gene_indices_2
