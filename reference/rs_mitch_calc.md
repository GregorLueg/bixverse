# Calculate mitch enrichment leveraging Rust under the hood

**\[experimental\]** Ranks the contrasts column-wise and runs the mitch
MANOVA/ANOVA tests per pathway.

## Usage

``` r
rs_mitch_calc(x, pathway_list, min_size)
```

## Arguments

- x:

  Numerical matrix. Each column represents one of the contrasts you wish
  to test for and the rows represent the gene statistics per contrast.
  Needs row names (the genes).

- pathway_list:

  Named list. Each element represents one of the pathways to test for.

- min_size:

  Integer. Minimum size of the gene set to be tested for.

## Value

A list with the following elements:

- pathway_names - The name of the pathway.

- pathway_sizes - The size of the pathway.

- manova_pval - The p-value of the MANOVA test.

- manova_fdr - The Benjamini-Hochberg adjusted `manova_pval`.

- anova_pvals - The p-values of the ANOVA test on top of the MANOVA
  results. Total length = `ncol(x)` \* number of pathways,
  pathway-major.

- scores - The scores for each pathway set, contrast. Same length and
  layout as `anova_pvals`.

- s_dist - Calculated distances from the hypotenuse.

- sd - SDs of the scores.
