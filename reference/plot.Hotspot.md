# Plot the Hotspot Z-score matrix

Produces a heatmap of the pairwise gene-gene Z-score matrix. If module
membership has been computed via
[`generate_hotspot_membership()`](https://gregorlueg.github.io/bixverse/reference/generate_hotspot_membership.md),
genes are ordered by module and unassigned genes are excluded. If no
membership is available, all genes are shown.

## Usage

``` r
# S3 method for class 'Hotspot'
plot(x, max_genes = 500L, seed = 42L, ...)
```

## Arguments

- x:

  A `Hotspot` object.

- max_genes:

  Integer. Maximum number of genes to plot. If the number of genes
  exceeds this, a random subsample is drawn (stratified by module if
  membership is available). Set to `NULL` to disable. Default 500.

- seed:

  Integer. Seed for reproducible subsampling.

- ...:

  Further arguments (currently unused).

## Value

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.
