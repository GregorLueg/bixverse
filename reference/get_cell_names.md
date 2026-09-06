# Get the cell names

Returns the cell names (usually barcodes).

## Usage

``` r
get_cell_names(x, filtered = FALSE)

# S3 method for class 'ScMap'
get_cell_names(x, filtered = FALSE)

## S7 method for class <bixverse::SingleCells>
get_cell_names(x, filtered = FALSE)

## S7 method for class <bixverse::SingleCellsSubset>
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

## Value

The cell names (barcodes)

## Examples

``` r
# barcodes of the cells that passed quality control
sc <- demo_single_cells(prepped = FALSE)
head(get_cell_names(sc, filtered = TRUE), 3)
#> [1] "cell_001" "cell_002" "cell_003"

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
