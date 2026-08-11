# Extract normalised gene expression for plotting

Extracts dense normalised (log1p) expression values for a set of genes,
optionally with additional observation metadata columns.

## Usage

``` r
extract_gene_expression(
  object,
  features,
  obs_cols = NULL,
  scale = FALSE,
  clip = NULL,
  modality = c("rna", "adt"),
  layer = c("norm", "magic")
)
```

## Arguments

- object:

  A single cell class.

- features:

  Character vector. Gene IDs to extract.

- obs_cols:

  Optional character vector. Column names from the obs table to include.

- scale:

  Boolean. Whether to z-score the expression values.

- clip:

  Optional numeric. If `scale = TRUE`, clip z-scores to `[-clip, clip]`.

- modality:

  String. One of `c("rna", "adt")`. ADT is only available for
  `SingleCellsMultiModal`.

- layer:

  String. One of `c("norm", "magic")`. With `"magic"` the values come
  from the imputed layer
  [`run_magic_sc()`](https://gregorlueg.github.io/bixverse/reference/run_magic_sc.md)
  wrote, which only holds the genes it was asked for. Imputation
  inflates gene-gene correlation, so this is for looking at things, not
  for measuring them. Note that
  [`extract_dot_plot_data()`](https://gregorlueg.github.io/bixverse/reference/extract_dot_plot_data.md)
  deliberately has no such argument: group means of imputed values are
  exactly the quantity MAGIC manufactures.

## Value

A data.table with a `cell_id` column, one column per gene, and any
requested obs columns.
