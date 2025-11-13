# Helper function to generate fgsea simple-based permutations

Helper function to generate fgsea simple-based permutations

## Usage

``` r
rs_calc_gsea_stat_cumulative_batch(
  stats,
  pathway_scores,
  pathway_sizes,
  iters,
  gsea_param,
  return_add_stats,
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

- gsea_param:

  Float. The Gene Set Enrichment parameter.

- return_add_stats:

  Boolean. Returns additional statistics necessary for the multi-level
  calculations.

- seed:

  Integer. For reproducibility purposes

## Value

List with the following elements

- es Enrichment scores for the gene sets

- nes Normalised enrichment scores for the gene sets

- pvals The calculated p-values.

- n_more_extreme Number of times the enrichment score was bigger or
  smaller than the permutation (pending sign).

- size Pathway size.

If `return_add_stats` is set to true, there is additional elements in
the list:

- le_zero Number of times the permutation was less than zero.

- ge_zero Number of times the permutation was greater than zero.
