# SCENIC: Select the Top TF \<\> Gene pairs

SCENIC: Select the Top TF \<\> Gene pairs

## Usage

``` r
rs_top_k_targets(matrix, k, margin, min_value)
```

## Arguments

- matrix:

  Numeric matrix with genes x TF importance values

- k:

  Integer. Number of top genes / TFs to extract.

- margin:

  If set to 1, the top k TFs per gene are used. If set to 2, the top k
  genes per TF are used. Both versions were used in the original paper.

- min_value:

  Float. An

## Value

A list with three vectors: tf, gene, importance
