# Get the cells to keep

Returns the indices of the cells that survived quality control. These
are stored 0-indexed for Rust, so add one before using them in R.

## Usage

``` r
get_cells_to_keep(x)

# S3 method for class 'ScMap'
get_cells_to_keep(x)

## S7 method for class <bixverse::SingleCells>
get_cells_to_keep(x)

## S7 method for class <bixverse::SingleCellsSubset>
get_cells_to_keep(x)
```

## Arguments

- x:

  An object from which to get the cells to keep from. These are
  0-indexed.

## Value

Integer vector with 0-indices of the cells to keep.

## Examples

``` r
# 0-based for Rust, so add one before indexing in R
sc <- demo_single_cells(prepped = FALSE)
head(get_cells_to_keep(sc) + 1, 3)
#> [1] 1 2 3

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
