# Calculate the pairwise gene-correlation for meta cells

**\[experimental\]** Correlates `gene_indices_1[i]` against
`gene_indices_2[i]` over the meta cells, in memory.

## Usage

``` r
rs_pairwise_gene_cors_mc(
  sparse_data,
  gene_indices_1,
  gene_indices_2,
  spearman,
  verbose
)
```

## Arguments

- sparse_data:

  A named list that needs to have `data`, `indptr`, `indices`, `nrow`,
  `ncol` and `cs_type`. Shape is (metacells, genes), holding the
  normalised counts.

- gene_indices_1:

  Integer. The gene indices for the first set of genes. Must be
  0-indexed!

- gene_indices_2:

  Integer. The gene indices for the second set of genes, same length as
  `gene_indices_1`. Must be 0-indexed!

- spearman:

  Boolean. Shall the Spearman correlation be calculated.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

Numeric vector with one correlation per pair of `gene_indices_1` and
`gene_indices_2`.
