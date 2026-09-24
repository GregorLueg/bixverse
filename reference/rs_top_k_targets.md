# SCENIC: Select the Top TF \<\> Gene pairs

**\[experimental\]**

## Usage

``` r
rs_top_k_targets(matrix, k, margin, min_value)
```

## Arguments

- matrix:

  Numeric matrix with genes (rows) x TFs (columns) importance values.
  Must carry row and column names.

- k:

  Integer. Number of top genes / TFs to extract.

- margin:

  Integer. If set to 1, the top k TFs per gene are used. If set to 2,
  the top k genes per TF are used. Both versions were used in the
  original paper. Any other value errors.

- min_value:

  Optional float. Pairs with an importance below this are never
  selected.

## Value

A list with three vectors: `tf`, `gene`, `importance`.
