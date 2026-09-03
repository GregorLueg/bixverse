# Plot the Hotspot Z-score matrix

Heatmap of the pairwise gene-gene Z-scores. With module membership from
[`generate_hotspot_membership()`](https://gregorlueg.github.io/bixverse/reference/generate_hotspot_membership.md)
the plot keeps the `top_k` largest modules, orders genes within each
module by hierarchical clustering and separates and labels the module
blocks, so the block structure is actually readable. Without membership
every gene is shown in the order it comes in.

## Usage

``` r
# S3 method for class 'Hotspot'
plot(x, top_k = 5L, max_genes = 500L, seed = 42L, ...)
```

## Arguments

- x:

  A `Hotspot` object.

- top_k:

  Integer. Number of modules to keep, ranked by gene count. Set to
  `NULL` to keep all of them. Ignored when no membership has been
  computed. Defaults to `5L`.

- max_genes:

  Integer. Maximum number of genes to plot. Above this a subsample is
  drawn, allocated across the kept modules in proportion to their size.
  Set to `NULL` to disable. Defaults to `500L`.

- seed:

  Integer. Seed for reproducible subsampling.

- ...:

  Further arguments (currently unused).

## Value

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.
