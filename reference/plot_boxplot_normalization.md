# Helper plot function for boxplot of normalized data

Helper plot function for boxplot of normalized data

## Usage

``` r
plot_boxplot_normalization(samples, voom_object, group_col)
```

## Arguments

- samples:

  data.table with sample information with perc_detected_genes and a
  column specifying the cohort.

- voom_object:

  `EList`. Voom object with normalised counts.

- group_col:

  String. The grouping column.

## Value

ggplot object, i.e., box plot with expression per sample.
