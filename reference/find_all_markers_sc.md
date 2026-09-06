# Find all markers

This function can be used to run differential gene expression for every
group of an unsupervised clustering method for example. You specify a
column and the function will start calculating differential gene
expression of the first cluster vs. everything else, second cluster vs.
everything else, etc. The function will automatically downsample
everything else to a random set of 100,000 cells if it should exceed
that. This automatic downsampling can be turned off however.

## Usage

``` r
find_all_markers_sc(
  object,
  column_of_interest,
  method = "wilcox",
  alternative = c("greater", "less", "twosided"),
  min_prop = 0.05,
  downsampling = TRUE,
  seed = 42L,
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

  String. Which method to use for the calculations of the DGE. At the
  moment the only option is `"wilcox"`, but the parameter is reserved
  for future features.

- alternative:

  String. Test alternative. One of `c("twosided", "greater", "less")`.
  This function will default to `"greater"`, i.e., genes upregulated in
  the group.

- min_prop:

  Numeric. The minimum proportion of cells that need to express the gene
  to be tested in any of the two groups.

- downsampling:

  Boolean. If the other group exceeds 100,000 cells, a random subsample
  of 100,000 cells will be used.

- seed:

  Integer. Seed that is used for the downsampling.

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
res <- find_all_markers_sc(
  sc,
  column_of_interest = "leiden_clustering",
  .verbose = FALSE
)
head(res)
#>      grp gene_id      lfc     prop1     prop2 z_scores     p_values
#>    <int>  <char>    <num>     <num>     <num>    <num>        <num>
#> 1:     0 gene_01 3.211139 0.9881657 0.7280967 16.50891 1.582561e-61
#> 2:     0 gene_02 3.121258 0.9881657 0.7129909 16.00064 6.323246e-58
#> 3:     0 gene_03 3.216300 0.9881657 0.7311178 16.73337 3.744067e-63
#> 4:     0 gene_04 3.098137 0.9644970 0.6706949 15.48991 2.029132e-54
#> 5:     0 gene_05 2.711545 0.8284023 0.4290030 12.43032 8.944763e-36
#> 6:     0 gene_06 3.202754 0.9822485 0.6767372 16.24377 1.236194e-59
#>             fdr
#>           <num>
#> 1: 1.978202e-60
#> 2: 5.269372e-57
#> 3: 9.360168e-62
#> 4: 1.449380e-53
#> 5: 4.969313e-35
#> 6: 1.236194e-58

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
