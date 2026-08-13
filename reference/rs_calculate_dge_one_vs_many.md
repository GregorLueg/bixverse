# Calculate one-vs-many AUROC DGEs for specific markers

**\[experimental\]** The function scores one reference group of cells
against each comparison group separately and summarises the results per
gene across all of the comparisons. This is the marker question: a gene
that is specific to the reference has to hold up against every rival,
which a single pooled test cannot answer because it is dominated by
whichever rival contributes the most cells. Genes are filtered once,
globally, so every comparison's FDR is calculated over the same gene
set.

## Usage

``` r
rs_calculate_dge_one_vs_many(
  f_path,
  cell_indices_ref,
  cell_indices_other,
  min_prop,
  alternative,
  verbose
)
```

## Arguments

- f_path:

  String. Path to the `counts_cells.bin` file.

- cell_indices_ref:

  Integer. Index positions (0-indexed) of the cells of the reference
  group.

- cell_indices_other:

  List. List of integer vectors, each containing the index positions
  (0-indexed) of the cells of one comparison group.

- min_prop:

  Minimum proportion of expression in at least one of the groups to be
  tested.

- alternative:

  String. One of `c("twosided", "greater", "less")`. Null hypothesis.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with the elements below. The per-comparison elements are
flattened comparison-major, i.e., all genes of the first comparison,
then all genes of the second, and so on.

- comparison - Index (0-indexed) of the comparison group the
  per-comparison statistics belong to.

- auroc - AUROC of the reference against the comparison group.

- lfc - Log fold change of the reference against the comparison group.

- prop_other - Proportion of cells expressing the gene in the comparison
  group.

- z_scores - Z-scores based on the Mann Whitney statistic.

- p_values - P-values of the Mann Whitney statistic.

- fdr - False discovery rate after BH adjustment, per comparison.

- prop_ref - Proportion of reference cells expressing the gene.

- median_auroc - Median AUROC across the comparisons.

- min_auroc - Worst AUROC across the comparisons.

- mean_auroc - Mean AUROC across the comparisons.

- max_auroc - Best AUROC across the comparisons.

- worst_comparison - Index (0-indexed) of the comparison group achieving
  `min_auroc`.

- min_rank - Best rank the gene achieves in any single comparison when
  the genes are ordered by descending AUROC.

- simes_p - Simes-combined p-value across the comparisons.

- simes_fdr - False discovery rate over `simes_p`.

- max_p - Largest p-value across the comparisons.

- max_p_fdr - False discovery rate over `max_p`.

- genes_to_keep - Boolean indicating which genes were tested.
