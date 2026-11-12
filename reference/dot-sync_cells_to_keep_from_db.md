# Take the cell filter from DuckDB when loading a saved object

DuckDB `obs.to_keep` is the authoritative cell filter: another front end
(or a hand edit of the DB) may have changed it after the memory file was
written. When the two disagree, the DB wins and the cached artefacts in
the memory file are flagged as stale.

## Usage

``` r
.sync_cells_to_keep_from_db(object, duckdb_con)
```

## Arguments

- object:

  `SingleCells` or `SingleCellsMultiModal` with the memory file already
  loaded.

- duckdb_con:

  The `SingleCellDuckDB` connection of the object.

## Value

The object with the cell filter set from DuckDB.
