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

  `single_cell_exp` class.

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

  Boolean. Controls verbosity of the function.

## Value

data.table with the DGE results from the test.
