# Cluster the genes by Z-score together

**\[experimental\]** Builds an average-linkage dendrogram over the pair
Z-scores and assigns modules as HotSpot's `compute_modules` does.

## Usage

``` r
rs_hotspot_cluster_genes(z_matrix, fdr_threshold, min_size)
```

## Arguments

- z_matrix:

  Numerical matrix. Symmetric gene x gene Z-scores with a zero diagonal,
  as returned by
  [`rs_hotspot_gene_cor()`](https://gregorlueg.github.io/bixverse/reference/rs_hotspot_gene_cor.md).
  Must be finite.

- fdr_threshold:

  Float. BH level at which a pair Z-score counts as significant.

- min_size:

  Integer. Minimum number of genes per module.

## Value

Numeric vector with one module label per gene, 0-indexed and numbered
densely. `NaN` indicates that the gene did not pass the thresholds and
has not been assigned.
