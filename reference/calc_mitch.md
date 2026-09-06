# Calculate a mitch gene set enrichments on contrast

Rust-based version of the mitch multi-contrast enrichment, see Kaspi and
Ziemann. Takes in a matrix representing the contrasts you wish to test
against.

## Usage

``` r
calc_mitch(contrast_mat, gene_set_list, min_size = 5L)
```

## Arguments

- contrast_mat:

  Numerical matrix. The rows represent the gene statistic per contrast
  and each column represents the contrast.

- gene_set_list:

  Named list. Contains the pathways you wish to test against.

- min_size:

  Integer. Minimum size of the gene set to be included.

## Value

A data.table with the Mitch enrichment results.

## References

Kaspi and Ziemann, Bmc Genomics, 2020

## Examples

``` r
# multi-contrast enrichment over two contrasts
set.seed(10L)
contrast_mat <- matrix(
  rnorm(300 * 2),
  ncol = 2,
  dimnames = list(sprintf("gene_%03i", 1:300), c("contrast_a", "contrast_b"))
)
contrast_mat[1:30, ] <- contrast_mat[1:30, ] + 1.5
gene_sets <- list(
  hit_set = sprintf("gene_%03i", 1:30),
  bg_set = sprintf("gene_%03i", 100:150)
)
res <- calc_mitch(contrast_mat, gene_sets)
res[, c("pathway_names", "manova_pval", "s_dist")]
#>    pathway_names manova_pval     s_dist
#>           <char>       <num>      <num>
#> 1:       hit_set   0.0000000 1.02835712
#> 2:        bg_set   0.8558181 0.04953152
```
