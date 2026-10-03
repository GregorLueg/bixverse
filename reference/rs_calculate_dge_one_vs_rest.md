# Calculate one-vs-rest Mann Whitney DGEs for every cell group

**\[experimental\]** Tests every cell group against all other grouped
cells pooled, in a single pass over the gene-based file. Cells in no
group are ignored. Each group filters genes on its own: a gene is tested
for a group if it clears `min_prop` in the group or in the rest, and the
FDR is calculated over that group's tested genes.

## Usage

``` r
rs_calculate_dge_one_vs_rest(
  f_path,
  cell_groups,
  min_prop,
  alternative,
  verbose
)
```

## Arguments

- f_path:

  String. Path to the `counts_genes.bin` file.

- cell_groups:

  List. List of integer vectors, each containing the index positions
  (0-indexed) of the cells of one group.

- min_prop:

  Numeric. Minimum proportion of expression in the group or in the rest
  to be tested.

- alternative:

  String. One of `c("twosided", "greater", "less")`. Null hypothesis.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with the elements below, one entry per tested gene and group,
flattened group-major.

- group - Index (0-indexed) of the group in `cell_groups`.

- gene_idx - Index (0-indexed) of the gene.

- lfc - Log fold change of the group against the rest.

- prop1 - Proportion of cells expressing the gene in the group.

- prop2 - Proportion of cells expressing the gene in the rest.

- auroc - AUROC of the group against the rest.

- z_scores - Z-scores based on the Mann Whitney statistic.

- p_values - P-values of the Mann Whitney statistic.

- fdr - False discovery rate after BH adjustment, per group.
