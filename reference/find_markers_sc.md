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
  method = NULL,
  alternative = c("twosided", "greater", "less"),
  min_prop = 0.05,
  modality = c("rna", "atac"),
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

  Optional string. One of `c("wilcox", "lr")`: the Wilcoxon rank sum
  test, or Seurat's logistic regression test (`LRDETest()`) with the
  library size as the latent variable. `NULL` picks `"lr"` on ATAC data
  and `"wilcox"` otherwise.

- alternative:

  String. Test alternative. One of `c("twosided", "greater", "less")`.
  Function will default to `"twosided"`.

- min_prop:

  Numeric. The minimum proportion of cells that need to express the gene
  to be tested in any of the two groups.

- modality:

  String. One of `c("rna", "atac")`. `"atac"` tests the features of the
  ATAC store of a `SingleCellsMultiModal`; on a `SingleCellsAtac` both
  name its own store. Results are then keyed by `feature_id`.

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
#> 1: gene_01 3.339092 0.9880239 0.6946108 14.38336 6.581099e-47 7.105061e-46
#> 2: gene_02 3.292497 0.9880239 0.6886228 14.16017 1.615884e-45 1.009927e-44
#> 3: gene_03 3.461720 1.0000000 0.6946108 14.95812 1.378589e-50 3.446472e-49
#> 4: gene_04 2.915461 0.9580838 0.6886228 12.95331 2.250255e-38 8.654827e-38
#> 5: gene_05 2.709359 0.8383234 0.4431138 10.69039 1.128701e-26 3.319709e-26
#> 6: gene_06 3.113653 0.9760479 0.7005988 14.03240 9.874756e-45 5.482066e-44

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
