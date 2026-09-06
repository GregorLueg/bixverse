# Set cell mapping

Set a cell mapping for a given object. This is used for the single
cell-related classes with streaming from disk.

## Usage

``` r
set_cell_mapping(x, cell_map)

# S3 method for class 'ScMap'
set_cell_mapping(x, cell_map)

## S7 method for class <bixverse::SingleCells>
set_cell_mapping(x, cell_map)
```

## Arguments

- x:

  An object to set cell mapping for

- cell_map:

  Named integer indicating indices and names of the cells

## Examples

``` r
# the mapping is normally written during ingestion
sc <- demo_single_cells(prepped = FALSE)
cells <- get_cell_names(sc)
sc <- set_cell_mapping(sc, stats::setNames(seq_along(cells), cells))
head(get_cell_names(sc), 3)
#> [1] "cell_001" "cell_002" "cell_003"

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
