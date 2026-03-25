# Get the cell names

Get the cell names

Get the cell names from a `SingleCells`.

## Usage

``` r
get_cell_names(x, filtered = FALSE)

# S3 method for class 'ScMap'
get_cell_names(x, filtered = FALSE)
```

## Arguments

- x:

  An object to get the cell names from.

- filtered:

  Boolean. Shall, if found only the cell names of the `cells_to_keep` be
  returned (see
  [`set_cells_to_keep()`](https://gregorlueg.github.io/bixverse/reference/set_cells_to_keep.md).
  Defaults to `FALSE`
