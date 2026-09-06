# Helper plot function for boxplot of normalised data

Helper plot function for boxplot of normalised data

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

## Examples

``` r
# normalised expression per sample, coloured by cohort
syn <- synthetic_bulk_cor_matrix()
samples <- data.table::data.table(
  sample_id = colnames(syn$counts),
  cohort = rep(c("case", "control"), each = 50)
)
dge_list <- edgeR::normLibSizes(edgeR::DGEList(counts = syn$counts))
voom_obj <- limma::voom(
  dge_list,
  stats::model.matrix(~ samples$cohort)
)
plot_boxplot_normalization(samples, voom_obj, group_col = "cohort")
```
