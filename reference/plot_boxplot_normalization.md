# Helper plot function for boxplot of normalised data

Helper plot function for boxplot of normalised data

## Usage

``` r
plot_boxplot_normalization(samples, norm_counts, group_col)
```

## Arguments

- samples:

  data.table with sample information with perc_detected_genes and a
  column specifying the cohort.

- norm_counts:

  Numeric matrix. The normalised log2 expression, genes x samples, with
  the samples in the same order as the rows of `samples`.

- group_col:

  String. The grouping column.

## Value

ggplot object, i.e., box plot with expression per sample.

## Examples

``` r
# normalised expression per sample, coloured by cohort
syn <- synthetic_bulk_cor_matrix()
samples <- data.table::data.table(
  sample_id = colnames(syn$counts),
  cohort = rep(c("case", "control"), each = 50)
)
norm_counts <- rs_cpm(syn$counts, lib_size = NULL, log = TRUE,
  prior_count = 0.5)
plot_boxplot_normalization(samples, norm_counts, group_col = "cohort")
```
