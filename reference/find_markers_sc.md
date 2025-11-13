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

  `single_cell_exp` class.

- cells_1:

  String. The names of the cells in group 1. Need to be part of the cell
  names in the object, see [`get_cell_names()`](get_cell_names.md).

- cells_2:

  String. The names of the cells in group 2. Need to be part of the cell
  names in the object, see [`get_cell_names()`](get_cell_names.md).

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

  Boolean. Controls verbosity of the function.

## Value

data.table with the DGE results from the test.
