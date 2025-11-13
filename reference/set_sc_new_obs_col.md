# Add a new column to the obs table

Add a new column to the obs table

## Usage

``` r
set_sc_new_obs_col(object, col_name, new_data)
```

## Arguments

- object:

  [`bixverse::single_cell_exp`](single_cell_exp.md) class.

- col_name:

  String. The name of the column to add.

- new_data:

  Atomic vector. The data to add to the column. Needs to be of same
  length as [`get_cells_to_keep()`](get_cells_to_keep.md) and have the
  same order.

## Value

The class with updated obs table in the DuckDB
