# Set cells to keep

Set the cells to keep. This is used for the single cell-related classes
with streaming from disk and tells subsequent (Rust) methods which cells
to include.

## Usage

``` r
set_cells_to_keep(x, cells_to_keep)

# S3 method for class 'ScMap'
set_cells_to_keep(x, cells_to_keep)

## S7 method for class <bixverse::SingleCells>
set_cells_to_keep(x, cells_to_keep)
```

## Arguments

- x:

  An object to set cells to keep for

- cells_to_keep:

  String or integer. The names or indices of the cells to keep in
  downstream analysis.

## Examples

``` r
# restrict everything downstream to the first 100 cells
sc <- demo_single_cells(prepped = FALSE)
sc <- set_cells_to_keep(sc, get_cell_names(sc)[1:100])
length(get_cells_to_keep(sc))
#> [1] 100

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
