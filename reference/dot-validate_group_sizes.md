# Validate that all groups meet minimum cell count requirements

Stops if any group has fewer than 50 cells. Warns if any group has fewer
than 500 cells.

## Usage

``` r
.validate_group_sizes(groups)
```

## Arguments

- groups:

  Named list of cell index vectors as returned by
  [`.split_cells_by_group()`](https://gregorlueg.github.io/bixverse/reference/dot-split_cells_by_group.md).

## Value

Invisibly `NULL`.
