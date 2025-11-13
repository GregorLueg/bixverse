# Helper plot function of distribution of genes by samples

Helper plot function of distribution of genes by samples

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
