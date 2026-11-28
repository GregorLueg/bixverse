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
  expr_modality = c("rna", "adt", "atac"),
  embd_modality = c("rna", "adt", "atac", "wnn"),
  layer = c("norm", "magic"),
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
  `c("rna", "adt", "atac")`.

- embd_modality:

  String. Modality the embedding is pulled from. One of
  `c("rna", "adt", "atac", "wnn")`. Use `"wnn"` for WNN-derived
  embeddings.

- layer:

  String. One of `c("norm", "magic")`, forwarded to
  [`extract_gene_expression()`](https://gregorlueg.github.io/bixverse/reference/extract_gene_expression.md).
  Use `"magic"` to colour the embedding by the imputed layer
  [`run_magic_sc()`](https://gregorlueg.github.io/bixverse/reference/run_magic_sc.md)
  wrote.

- ...:

  Additional arguments forwarded to
  [`extract_embedding_data()`](https://gregorlueg.github.io/bixverse/reference/extract_embedding_data.md)
  and onward to
  [`get_embedding()`](https://gregorlueg.github.io/bixverse/reference/get_embedding.md).
  Do not pass `modality` here; the embedding modality is set via
  `embd_modality` and passing it again will error.

## Value

A long data.table with `cell_id`, `dim_*`, `gene` and `expression`.

## Examples

``` r
# two genes melted onto the PCA coordinates
sc <- demo_single_cells()
dt <- extract_feature_plot_data(
  sc,
  features = get_gene_names(sc)[1:2],
  embedding = "pca"
)
head(dt[, c("cell_id", "dim_1", "dim_2", "gene", "expression")])
#>     cell_id      dim_1      dim_2    gene expression
#>      <char>      <num>      <num>  <fctr>      <num>
#> 1: cell_001  0.3693915  3.0684264 gene_01   6.223760
#> 2: cell_002 -2.2334790 -2.5043652 gene_01   4.111858
#> 3: cell_003  2.4827571  0.5070686 gene_01   4.299178
#> 4: cell_004 -1.3581411  2.6233981 gene_01   5.756331
#> 5: cell_005 -2.1339433 -0.2906340 gene_01   0.000000
#> 6: cell_006  2.6228275  0.5391235 gene_01   0.000000

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
