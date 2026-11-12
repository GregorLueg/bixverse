# Write the object's cell filter to DuckDB before saving

R copies of a `SingleCells` object share one DuckDB, so a
[`set_cells_to_keep()`](https://gregorlueg.github.io/bixverse/reference/set_cells_to_keep.md)
on one copy rewrites `obs.to_keep` under all others. Saving re-asserts
the saved object's own filter, so the memory file and the DB agree on
disk.

## Usage

``` r
.sync_cells_to_keep_to_db(object)
```

## Arguments

- object:

  `SingleCells` or `SingleCellsMultiModal` class.

## Value

Invisible `NULL`.
