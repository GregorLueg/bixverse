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
