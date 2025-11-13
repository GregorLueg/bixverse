# Calculate mitch enrichment leveraging Rust under the hood

Calculate mitch enrichment leveraging Rust under the hood

## Usage

``` r
rs_mitch_calc(x, pathway_list, min_size)
```

## Arguments

- x:

  Numerical matrix. Each column represents on the contrasts you wish to
  test for and the rows represent the gene statistics per contrast.

- pathway_list:

  Named list. Each element represents one of the pathways to test for.

- min_size:

  Integer. Minimum size of gene the gene set to be tested for.

## Value

A list with the following elements:

- pathway_names - The name of the pathway.

- pathway_sizes The size of the pathway.

- manova_pvals - The p-value of the MANOVA test.

- anova_pvals The p-values of the ANOVA test on top of the MANOVA
  results. Total length = `ncol(x)` \* number of pathways.

- scores - The scores for each pathway set, contrast. Same length as
  `anova_pvals`.

- s_dist - Calculated distances from the hypotenuse.

- sd - SDs of the scores.
