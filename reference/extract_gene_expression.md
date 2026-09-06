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

## Examples

``` r
# normalised expression of three genes with a cell annotation
sc <- demo_single_cells()
dt <- extract_gene_expression(
  sc,
  features = get_gene_names(sc)[1:3],
  obs_cols = "cell_grp"
)
head(dt, 3)
#>     cell_id  gene_01  gene_02  gene_03    cell_grp
#>      <char>    <num>    <num>    <num>      <char>
#> 1: cell_001 6.222656 6.761719 5.531250 cell_type_1
#> 2: cell_002 4.113281 4.511719 4.113281 cell_type_2
#> 3: cell_003 4.300781 4.300781 3.226562 cell_type_3

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
