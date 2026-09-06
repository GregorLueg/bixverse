# Rename columns in the obs or var table

Renames the columns in the obs or var table of single cell-related
classes.

## Usage

``` r
setnames_sc(object, table = c("obs", "var"), old, new)
```

## Arguments

- object:

  `SingleCells`, `MetaCells` (or potentially other) class.

- table:

  String. One of `c("obs", "var")`. In which of the tables to rename the
  columns.

- old:

  Character vector. The old column names.

- new:

  Character vector. The new column names.

## Value

Invisible self

## Examples

``` r
# rename a column in the obs table
sc <- demo_single_cells(prepped = FALSE)
sc <- setnames_sc(sc, table = "obs", old = "cell_grp", new = "cell_type")
head(get_sc_obs(sc)$cell_type, 3)
#> [1] "cell_type_1" "cell_type_2" "cell_type_3"

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
