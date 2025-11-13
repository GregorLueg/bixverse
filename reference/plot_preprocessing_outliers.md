# Helper plot function for identification of outliers

Helper plot function for identification of outliers

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
