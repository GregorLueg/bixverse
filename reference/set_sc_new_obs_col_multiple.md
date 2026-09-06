# Add multiple new columns to the obs table

Add multiple new columns to the obs table

## Usage

``` r
set_sc_new_obs_col_multiple(object, new_data)
```

## Arguments

- object:

  `SingleCells` class.

- new_data:

  Named list. The names will be the column names and the elements will
  be added to the obs table. Needs to be of same length as
  [`get_cells_to_keep()`](https://gregorlueg.github.io/bixverse/reference/get_cells_to_keep.md)
  and have the same order!

## Value

The class with updated obs table in the DuckDB

## Examples

``` r
# several columns in one go
sc <- demo_single_cells(prepped = FALSE)
n_cells <- dim(sc)[1]
sc <- set_sc_new_obs_col_multiple(
  sc,
  new_data = list(
    arm = rep(c("ctrl", "treated"), length.out = n_cells),
    donor = rep(sprintf("donor_%i", 1:4), length.out = n_cells)
  )
)
head(get_sc_obs(sc)[, c("arm", "donor")], 3)
#>        arm   donor
#>     <char>  <char>
#> 1:    ctrl donor_1
#> 2: treated donor_2
#> 3:    ctrl donor_3

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
