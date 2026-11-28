# Find all markers

This function can be used to run differential gene expression for every
group of an unsupervised clustering method for example. You specify a
column and the function calculates differential gene expression of the
first cluster vs. everything else, second cluster vs. everything else,
etc. All groups are tested in a single pass over the genes, so memory
stays bounded however many cells there are.

## Usage

``` r
find_all_markers_sc(
  object,
  column_of_interest,
  method = NULL,
  alternative = c("greater", "less", "twosided"),
  min_prop = 0.05,
  modality = c("rna", "atac"),
  .verbose = TRUE
)
```

## Arguments

- object:

  `SingleCells` or `SingleCellsSubset` class.

- column_of_interest:

  String. The column you wish to use to identify the markers between all
  combination. Needs to be in the obs table

- method:

  Optional string. One of `c("wilcox", "lr")`: the Wilcoxon rank sum
  test, or Seurat's logistic regression test (`LRDETest()`) with the
  library size as the latent variable. `NULL` picks `"lr"` on ATAC data
  and `"wilcox"` otherwise.

- alternative:

  String. Test alternative. One of `c("twosided", "greater", "less")`.
  This function will default to `"greater"`, i.e., genes upregulated in
  the group.

- min_prop:

  Numeric. The minimum proportion of cells that need to express the gene
  to be tested in the group or in the rest.

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
# one versus rest across the Leiden clusters
sc <- demo_single_cells()
sc <- find_clusters_sc(sc, res = 1.0)
#> Running leiden clustering at resolution 1.
res <- find_all_markers_sc(
  sc,
  column_of_interest = "leiden_clustering",
  .verbose = FALSE
)
head(res)
#>      grp gene_id      lfc     prop1     prop2 z_scores     p_values
#>    <int>  <char>    <num>     <num>     <num>    <num>        <num>
#> 1:     0 gene_01 3.210998 0.9881657 0.7280967 16.50923 1.574405e-61
#> 2:     0 gene_02 3.121287 0.9881657 0.7129909 16.00128 6.258756e-58
#> 3:     0 gene_03 3.216173 0.9881657 0.7311178 16.73269 3.786896e-63
#> 4:     0 gene_04 3.098070 0.9644970 0.6706949 15.48957 2.039930e-54
#> 5:     0 gene_05 2.711567 0.8284023 0.4290030 12.43066 8.906978e-36
#> 6:     0 gene_06 3.202758 0.9822485 0.6767372 16.24343 1.243042e-59
#>             fdr
#>           <num>
#> 1: 1.968006e-60
#> 2: 5.215630e-57
#> 3: 9.467241e-62
#> 4: 1.457093e-53
#> 5: 4.948321e-35
#> 6: 1.243042e-58

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
