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
  modality = c("rna", "adt"),
  layer = c("norm", "magic")
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

- layer:

  String. One of `c("norm", "magic")`, forwarded to
  [`extract_gene_expression()`](https://gregorlueg.github.io/bixverse/reference/extract_gene_expression.md).

## Value

A long data.table with `cell_id`, `group`, `gene` and `expression`.
`gene` is an ordered factor following `features`.

## Examples

``` r
# long format expression grouped by cell type
sc <- demo_single_cells(prepped = FALSE)
head(extract_gene_violin_data(
  sc,
  features = get_gene_names(sc)[1:2],
  grouping_variable = "cell_grp"
))
#>     cell_id       group    gene expression
#>      <char>      <fctr>  <fctr>      <num>
#> 1: cell_001 cell_type_1 gene_01   6.222656
#> 2: cell_002 cell_type_2 gene_01   4.113281
#> 3: cell_003 cell_type_3 gene_01   4.300781
#> 4: cell_004 cell_type_1 gene_01   5.757812
#> 5: cell_005 cell_type_2 gene_01   0.000000
#> 6: cell_006 cell_type_3 gene_01   0.000000

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
