# Load an existing SingleCells from disk

Helper function that can load the parameters to access the on-disk
stored data into the class.

## Usage

``` r
load_existing(object, .verbose = TRUE)
```

## Arguments

- object:

  `SingleCells`, `MetaCells` or `SingleCellsMultiModal` class.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

The object with added information on the data on disk.

## Examples

``` r
# a fresh handle over a directory written earlier
sc <- demo_single_cells(prepped = FALSE)
dir <- sc@dir_data
sc <- load_existing(SingleCells(dir_data = dir), .verbose = FALSE)
dim(sc)
#> [1] 500  50

unlink(dir, recursive = TRUE, force = TRUE)
```
