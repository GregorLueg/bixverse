# Extract grouped gene statistics for dot plots

Extracts per-group mean expression and percentage of expressing cells
for a set of genes. Returns a long-format data.table suitable for dot
plots.

## Usage

``` r
extract_dot_plot_data(
  object,
  features,
  grouping_variable,
  scale_exp = TRUE,
  modality = c("rna", "adt")
)
```

## Arguments

- object:

  A single cell class.

- features:

  Character vector. Gene IDs to extract.

- grouping_variable:

  String. Column name in the obs table to group by.

- scale_exp:

  Boolean. Whether to min-max scale mean expression per gene.

- modality:

  String. One of `c("rna", "adt")`. ADT is only available for
  `SingleCellsMultiModal`.

## Value

A data.table with columns: gene, group, mean_exp, scaled_exp, pct_exp.
