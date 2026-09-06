# Helper plot function for identification of outliers

Beeswarm plot of the percentage of detected genes per sample, with the
outlier cutoffs drawn in. Used in the bulk DGE pre-processing report.

## Usage

``` r
plot_preprocessing_outliers(samples, group_col, min_perc, max_perc)
```

## Arguments

- samples:

  data.table with sample information with perc_detected_genes and a
  column specifying the cohort.

- group_col:

  String specifying the column with cohort information

- min_perc:

  Numeric. Lower cutoff to identify outliers.

- max_perc:

  Numeric. Upper cutoff to identify outliers

## Value

ggplot object, i.e., beeswarm plot with outlier indication

## Examples

``` r
# percentage of detected genes with the outlier cutoffs drawn in
syn <- synthetic_bulk_cor_matrix()
samples <- data.table::data.table(
  cohort = rep(c("case", "control"), each = 50),
  perc_detected_genes = colMeans(syn$counts > 0) * 100
)
plot_preprocessing_outliers(
  samples,
  group_col = "cohort",
  min_perc = 60,
  max_perc = 95
)
```
