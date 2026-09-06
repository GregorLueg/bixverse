# Add a new column to the obs table

Add a new column to the obs table

## Usage

``` r
set_sc_new_obs_col(object, col_name, new_data)
```

## Arguments

- object:

  `SingleCells` class.

- col_name:

  String. The name of the column to add.

- new_data:

  Atomic vector. The data to add to the column. Needs to be of same
  length as
  [`get_cells_to_keep()`](https://gregorlueg.github.io/bixverse/reference/get_cells_to_keep.md)
  and have the same order!

## Value

The class with updated obs table in the DuckDB

## Examples

``` r
# one new column, in the order of the cells that passed quality control
sc <- demo_single_cells(prepped = FALSE)
sc <- set_sc_new_obs_col(
  sc,
  col_name = "arm",
  new_data = rep(c("ctrl", "treated"), length.out = dim(sc)[1])
)
table(get_sc_obs(sc)$arm)
#> 
#>    ctrl treated 
#>     250     250 

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
