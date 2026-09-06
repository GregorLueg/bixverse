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

A data.table with columns: gene, group, mean_exp, scaled_exp and
pct_exp.

## Examples

``` r
# mean expression and expressing fraction per cell group
sc <- demo_single_cells()
dt <- extract_dot_plot_data(
  sc,
  features = get_gene_names(sc)[1:5],
  grouping_variable = "cell_grp"
)
head(dt, 3)
#>       gene       group mean_exp  pct_exp scaled_exp
#>     <fctr>      <fctr>    <num>    <num>      <num>
#> 1: gene_01 cell_type_1 6.241650 98.80239 1.00000000
#> 2: gene_01 cell_type_2 2.902390 69.46108 0.00000000
#> 3: gene_01 cell_type_3 3.123859 76.50602 0.06632256

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
