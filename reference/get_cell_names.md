# Get the cell names

Get the cell names

Get the cell names from a `single_cell_exp`.

## Usage

``` r
get_cell_names(x, filtered = FALSE)

# S3 method for class 'sc_mapper'
get_cell_names(x, filtered = FALSE)
```

## Arguments

- x:

  An object to get the cell names from.

- filtered:

  Boolean. Shall, if found only the cell names of the `cells_to_keep` be
  returned (see [`set_cells_to_keep()`](set_cells_to_keep.md). Defaults
  to `FALSE`
