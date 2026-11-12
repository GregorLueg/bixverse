# Attach obs columns to extracted rows by cell id

Joins by id rather than by position, so the obs values cannot slide off
the extracted values whatever order either side comes in.

## Usage

``` r
.attach_obs_by_id(dt, object, obs_cols)
```

## Arguments

- dt:

  data.table. Extracted rows with a `cell_id` column.

- object:

  `SingleCells`, `SingleCellsSubset` or `MetaCells` class.

- obs_cols:

  Character vector. The obs columns to attach.

## Value

`dt` with the obs columns added by reference.
