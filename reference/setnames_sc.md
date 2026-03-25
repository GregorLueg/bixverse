# Rename columns in the obs or var table

Renames the columns in the obs or var table of single cell-related
classes.

## Usage

``` r
setnames_sc(object, table = c("obs", "var"), old, new)
```

## Arguments

- object:

  `SingleCells`, `MetaCells` class.

- table:

  String. One of `c("obs", "var")`. In which of the tables to rename the
  columns.

- old:

  Character vector. The old column names.

- new:

  Character vector. The new column names.

## Value

Invisible self
