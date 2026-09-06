# Calculate DGE between two cell groups

This function can be used to calculate differentially expressed genes
between two groups in the single cell data. At the moment, it has only
an implementation for the Wilcox-based rank statistic.

## Usage

``` r
find_markers_sc(
  object,
  cells_1,
  cells_2,
  method = c("wilcox"),
  alternative = c("twosided", "greater", "less"),
  min_prop = 0.05,
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells` or `SingleCellsSubset` class.

- cells_1:

  String. The names of the cells in group 1. Need to be part of the cell
  names in the object, see
  [`get_cell_names()`](https://gregorlueg.github.io/bixverse/reference/get_cell_names.md).

- cells_2:

  String. The names of the cells in group 2. Need to be part of the cell
  names in the object, see
  [`get_cell_names()`](https://gregorlueg.github.io/bixverse/reference/get_cell_names.md).

- method:

  String. Which method to use for the calculations of the DGE. At the
  moment the only option is `"wilcox"`, but the parameter is reserved
  for future features.

- alternative:

  String. Test alternative. One of `c("twosided", "greater", "less")`.
  Function will default to `"twosided"`.

- min_prop:

  Numeric. The minimum proportion of cells that need to express the gene
  to be tested in any of the two groups.

- .verbose:

  Boolean or integer. Controls verbosity and returns run times. `FALSE`
  -\> quiet, `TRUE` or `1L` -\> normal verbosity, `2L` -\> detailed
  verbosity.

## Value

data.table with the DGE results from the test.

## Examples

``` r
# Wilcoxon test between two of the planted cell types
sc <- demo_single_cells()
obs <- get_sc_obs(sc)
res <- find_markers_sc(
  sc,
  cells_1 = obs$cell_id[obs$cell_grp == "cell_type_1"],
  cells_2 = obs$cell_id[obs$cell_grp == "cell_type_2"],
  .verbose = FALSE
)
head(res)
#>    gene_id      lfc     prop1     prop2 z_scores     p_values          fdr
#>     <char>    <num>     <num>     <num>    <num>        <num>        <num>
#> 1: gene_01 3.339259 0.9880239 0.6946108 14.38283 6.632426e-47 7.043991e-46
#> 2: gene_02 3.292560 0.9880239 0.6886228 14.16134 1.589201e-45 9.932503e-45
#> 3: gene_03 3.461920 1.0000000 0.6946108 14.95929 1.354582e-50 3.386456e-49
#> 4: gene_04 2.915536 0.9580838 0.6886228 12.95389 2.233137e-38 8.588988e-38
#> 5: gene_05 2.709311 0.8383234 0.4431138 10.68923 1.142846e-26 3.361312e-26
#> 6: gene_06 3.113609 0.9760479 0.7005988 14.03355 9.715047e-45 5.397248e-44

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
