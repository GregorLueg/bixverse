# Extract per-cell expression mapped onto an embedding

Combines
[`extract_gene_expression()`](https://gregorlueg.github.io/bixverse/reference/extract_gene_expression.md)
with
[`extract_embedding_data()`](https://gregorlueg.github.io/bixverse/reference/extract_embedding_data.md)
and melts to long format, ready for faceted feature plots. The
expression source and the embedding source are chosen independently via
`expr_modality` and `embd_modality`, so you can colour an embedding from
one modality by expression from another (e.g. RNA expression on an
ADT-derived UMAP, or either modality on a WNN embedding). All sources
key on the same kept-cell barcodes, so the merge stays aligned
regardless of the chosen combination.

## Usage

``` r
extract_feature_plot_data(
  object,
  features,
  embedding,
  scale = FALSE,
  clip = NULL,
  obs_col = NULL,
  expr_modality = c("rna", "adt"),
  embd_modality = c("rna", "adt", "wnn"),
  ...
)
```

## Arguments

- object:

  A single cell class.

- features:

  Character vector. Gene/feature IDs to extract, taken from
  `expr_modality`.

- embedding:

  String. Name of the embedding.

- scale:

  Boolean. Whether to z-score the expression values.

- clip:

  Optional numeric. Clip z-scores if `scale = TRUE`.

- obs_col:

  Optional character vector. Obs columns to attach.

- expr_modality:

  String. Modality the expression is pulled from. One of
  `c("rna", "adt")`.

- embd_modality:

  String. Modality the embedding is pulled from. One of
  `c("rna", "adt", "wnn")`. Use `"wnn"` for WNN-derived embeddings.

- ...:

  Additional arguments forwarded to
  [`extract_embedding_data()`](https://gregorlueg.github.io/bixverse/reference/extract_embedding_data.md)
  and onward to
  [`get_embedding()`](https://gregorlueg.github.io/bixverse/reference/get_embedding.md).
  Do not pass `modality` here; the embedding modality is set via
  `embd_modality` and passing it again will error.

## Value

A long data.table with `cell_id`, `dim_*`, `gene` and `expression`.
