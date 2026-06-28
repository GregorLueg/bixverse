# Get the cell idx (R-based) and cell names

Returns the cell indices (R-based) and the cell names (usually barcodes)
from the object for further downstream usage.

## Usage

``` r
get_cell_info(x, filtered = TRUE)

# S3 method for class 'ScMap'
get_cell_info(x, filtered = TRUE)

## S7 method for class <bixverse::SingleCells>
get_cell_info(x, filtered = TRUE)
```

## Arguments

- x:

  An object to get the cell info from

- filtered:

  Boolean. If `TRUE`, only the cells to keep will be returned.

## Value

A named vector with elements -\> cell_idx, names -\> cell_names.
