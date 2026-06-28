# Extract per-cell expression grouped for violin plots

Combines
[`extract_gene_expression()`](https://gregorlueg.github.io/bixverse/reference/extract_gene_expression.md)
with a grouping obs column and melts to long format, ready for stacked
(one gene per row) violin plots.

## Usage

``` r
extract_gene_violin_data(
  object,
  features,
  grouping_variable,
  scale = FALSE,
  clip = NULL,
  modality = c("rna", "adt")
)
```

## Arguments

- object:

  A single cell class.

- features:

  Character vector. Gene IDs to extract.

- grouping_variable:

  String. Obs column to group by.

- scale:

  Boolean. Whether to z-score the expression values.

- clip:

  Optional numeric. Clip z-scores if `scale = TRUE`.

- modality:

  String. One of `c("rna", "adt")`.

## Value

A long data.table with `cell_id`, `group`, `gene` and `expression`.
`gene` is an ordered factor following `features`.
