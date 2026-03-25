# Getter the obs table

Getter the obs table

## Usage

``` r
get_sc_obs(object, indices = NULL, cols = NULL, filtered = FALSE)
```

## Arguments

- object:

  `SingleCells`, `MetaCells` class.

- indices:

  Optional integer vector. The integer positions of the cells to return.

- cols:

  Optional string vector. The columns from the obs table to return.

- filtered:

  Boolean. Whether to return all cells or filtered to `to_keep` cells.
  Not relevant for `MetaCells`.

## Value

The obs table
