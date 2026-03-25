# Add multiple new columns to the obs table

Add multiple new columns to the obs table

## Usage

``` r
set_sc_new_obs_col_multiple(object, new_data)
```

## Arguments

- object:

  [`bixverse::SingleCells`](https://gregorlueg.github.io/bixverse/reference/SingleCells.md)
  class.

- new_data:

  Named list. The names will be the column names and the elements will
  be added to the obs table. Needs to be of same length as
  [`get_cells_to_keep()`](https://gregorlueg.github.io/bixverse/reference/get_cells_to_keep.md)
  and have the same order!

## Value

The class with updated obs table in the DuckDB
