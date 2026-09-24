# Run a single hypergeometric test.

**\[experimental\]** Given a set of target genes, this is a Rust
implementation of an hypergeometric test testing for overenrichment of
the target genes in the gene sets.

## Usage

``` r
rs_hypergeom_test(
  target_genes,
  gene_sets,
  gene_universe,
  min_overlap,
  fdr_threshold
)
```

## Arguments

- target_genes:

  String vector. Represents the target gene set.

- gene_sets:

  List. Contains the strings that represent the gene sets to test
  against.

- gene_universe:

  String vector. The features representing the gene universe from which
  the target genes and gene sets are sampled from.

- min_overlap:

  Optional integer. Shall a filter be applied on the minimum of
  overlapping genes.

- fdr_threshold:

  Optional float. Shall a filter be applied for the maximum tolerated
  FDR.

## Value

A list containing (only for the gene sets passing the filters):

- pvals - The p-values from the hypergeometric test

- odds_ratios - The calculated odds ratios

- hits - The size of the overlap

- gene_set_lengths - The length of the gene sets.

- fdr - The FDR calculated across all gene sets (before filtering).

- to_keep - 1-based indices of the gene sets that passed the (optional)
  thresholds.
