# Helper plot function of distribution of genes by samples

Boxplot of the number of detected genes per sample, split by cohort.
Used in the bulk DGE pre-processing report to spot samples with poor
library complexity.

## Usage

``` r
plot_preprocessing_genes(samples, group_col)
```

## Arguments

- samples:

  data.table with sample information with nb_detected_genes and a column
  specifying the cohort.

- group_col:

  String specifying the column with cohort information

## Value

ggplot object, i.e. boxplot with number of genes by cohort

## Examples

``` r
# detected genes per sample, split by cohort
syn <- synthetic_bulk_cor_matrix()
samples <- data.table::data.table(
  cohort = rep(c("case", "control"), each = 50),
  nb_detected_genes = colSums(syn$counts > 0)
)
plot_preprocessing_genes(samples, group_col = "cohort")
```
