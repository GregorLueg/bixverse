# Split cell indices into per-group lists

Filters the obs table to `cells_to_use` and splits the resulting
(0-indexed) cell indices by the levels of `group_by`.

## Usage

``` r
.split_cells_by_group(object, group_by, cells_to_use)
```

## Arguments

- object:

  A `SingleCells` object.

- group_by:

  Character. Column name to group by.

- cells_to_use:

  Integer vector of 0-indexed cell indices.

## Value

Named list of 0-indexed integer vectors, one element per group.
