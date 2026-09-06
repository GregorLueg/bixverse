# Getter the obs table

Getter the obs table

## Usage

``` r
get_sc_obs(object, indices = NULL, cols = NULL, filtered = FALSE)
```

## Arguments

- object:

  `SingleCells`, `MetaCells`, `SingleCellsMultiModal` class.

- indices:

  Optional integer vector. The integer positions of the cells to return.

- cols:

  Optional string vector. The columns from the obs table to return.

- filtered:

  Boolean. Whether to return all cells or filtered to `to_keep` cells.
  Not relevant for `MetaCells`.

## Value

The obs table

## Examples

``` r
# the obs table, restricted to a few columns
sc <- demo_single_cells(prepped = FALSE)
head(get_sc_obs(sc, cols = c("cell_id", "cell_grp", "lib_size")), 3)
#>     cell_id    cell_grp lib_size
#>      <char>      <char>    <num>
#> 1: cell_001 cell_type_1      278
#> 2: cell_002 cell_type_2      333
#> 3: cell_003 cell_type_3      413

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
