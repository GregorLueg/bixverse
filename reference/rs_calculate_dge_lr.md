# Calculate DGEs between cells with a logistic regression test

**\[experimental\]** Seurat's `LRDETest()` with the library size as the
latent variable: a logistic regression of group membership on the gene
and the library size, against one on the library size alone, compared by
a likelihood ratio test. The default marker test for ATAC data.

## Usage

``` r
rs_calculate_dge_lr(
  f_path_cell,
  f_path_gene,
  cell_indices_1,
  cell_indices_2,
  min_prop,
  alternative,
  verbose
)
```

## Arguments

- f_path_cell:

  String. Path to the `counts_cells.bin` file, for the library sizes.

- f_path_gene:

  String. Path to the `counts_genes.bin` file.

- cell_indices_1:

  Integer. Index positions (0-indexed) of the cells of group 1.

- cell_indices_2:

  Integer. Index positions (0-indexed) of the cells of group 2.

- min_prop:

  Numeric. Minimum proportion of expression in at least one of the two
  groups to be tested.

- alternative:

  String. One of `c("twosided", "greater", "less")`. Null hypothesis.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with the following elements

- lfc - Log fold changes between the two groups.

- prop1 - Proportion of cells expressing the gene in group 1.

- prop2 - Proportion of cells expressing the gene in group 2.

- z_scores - Signed root of the likelihood ratio statistic.

- p_values - P-values of the likelihood ratio test.

- fdr - False discovery rate after BH adjustment

- genes_to_keep - Boolean indicating which genes were tested.

## References

Ntranos, et al., Nat Methods, 2019
