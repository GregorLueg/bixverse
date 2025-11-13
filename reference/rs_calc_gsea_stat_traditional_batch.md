# Helper function to generate traditional GSEA-based permutations

Helper function to generate traditional GSEA-based permutations

## Usage

``` r
rs_calc_gsea_stat_traditional_batch(
  stats,
  pathway_scores,
  pathway_sizes,
  iters,
  seed
)
```

## Arguments

- stats:

  Numeric vector. The gene level statistic. Needs to sorted in
  descending nature.

- pathway_scores:

  Numeric vector. The enrichment scores for the pathways

- pathway_sizes:

  Integer vector. The sizes of the pathways.

- iters:

  Integer. Number of permutations.

- seed:

  Integer For reproducibility purposes

## Value

List with the following elements

- es Enrichment scores for the gene sets

- nes Normalised enrichment scores for the gene sets

- pvals The calculated p-values.

- n_more_extreme Number of times the enrichment score was bigger or
  smaller than the permutation (pending sign).

- size Pathway size.
