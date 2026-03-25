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
  clip = NULL
)
```

## Arguments

- object:

  `SingleCells` class.

- features:

  Character vector. Gene IDs to extract.

- obs_cols:

  Optional character vector. Column names from the obs table to include.

- scale:

  Boolean. Whether to z-score the expression values.

- clip:

  Optional numeric. If `scale = TRUE`, clip z-scores to `[-clip, clip]`.

## Value

A data.table with a `cell_id` column, one column per gene, and any
requested obs columns.
