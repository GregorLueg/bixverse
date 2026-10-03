# Calculate one-vs-many AUROC DGEs for specific markers

**\[experimental\]** The function scores each reference group of cells
against every other group separately and summarises the results per gene
across all of the comparisons. This is the marker question: a gene that
is specific to the reference has to hold up against every rival, which a
single pooled test cannot answer because it is dominated by whichever
rival contributes the most cells. All reference groups come out of a
single pass over the gene-based file. Genes are filtered once, globally,
so every comparison's FDR is calculated over the same gene set.

## Usage

``` r
rs_calculate_dge_one_vs_many(
  f_path,
  cell_groups,
  references,
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

- references:

  Integer. Index positions (0-indexed) into `cell_groups` of the
  reference groups to report. The rivals of each reference are all other
  groups, in the order of `cell_groups`.

- min_prop:

  Numeric. Minimum proportion of expression in at least one of the
  groups to be tested.

- alternative:

  String. One of `c("twosided", "greater", "less")`. Null hypothesis.

- verbose:

  Integer. `0L` - quiet; `1L` - normal verbosity; `2L` - detailed
  verbosity.

## Value

A list with the elements below. The per-comparison elements are
flattened reference-major, then rival-major, then gene. The summary
elements are flattened reference-major, then gene.

- reference - Index (0-indexed) of the reference group, per comparison
  row.

- rival - Index (0-indexed) of the rival group, per comparison row.

- auroc - AUROC of the reference against the rival.

- lfc - Log fold change of the reference against the rival.

- prop_other - Proportion of cells expressing the gene in the rival.

- z_scores - Z-scores based on the Mann Whitney statistic.

- p_values - P-values of the Mann Whitney statistic.

- fdr - False discovery rate after BH adjustment, per comparison.

- summary_reference - Index (0-indexed) of the reference group, per
  summary row.

- prop_ref - Proportion of reference cells expressing the gene.

- median_auroc - Median AUROC across the rivals.

- min_auroc - Worst AUROC across the rivals.

- mean_auroc - Mean AUROC across the rivals.

- max_auroc - Best AUROC across the rivals.

- worst_rival - Index (0-indexed) of the rival group achieving
  `min_auroc`.

- min_rank - Best rank the gene achieves against any single rival when
  the genes are ordered by descending AUROC.

- simes_p - Simes-combined p-value across the rivals.

- simes_fdr - False discovery rate over `simes_p`.

- max_p - Largest p-value across the rivals.

- max_p_fdr - False discovery rate over `max_p`.

- genes_to_keep - Boolean indicating which genes were tested.
