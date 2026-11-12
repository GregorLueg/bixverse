# Load a plain-text metadata file into DuckDB

Reads a CSV / TSV (delimiter from the extension) into `table`, renaming
the first column to `id_col` and snake_casing the others when there is a
header, and numbering the kept rows in `idx_col`.

Row identity has to follow the file, and neither a window function over
`read_csv()` nor an `IN` list with one bound parameter per kept row is
safe for that at scale. The file is therefore first materialised into a
temporary table with insertion order preserved, whose `rowid` is then
the 0-based file row. The filter goes in as a registered data.frame and
a semi join, and the index is assigned ordered by `rowid`.

## Usage

``` r
.plain_text_to_duckdb(con, f_path, has_hdr, filter, table, idx_col, id_col)
```

## Arguments

- con:

  DuckDB connection.

- f_path:

  String. The file.

- has_hdr:

  Boolean. Whether the file has a header.

- filter:

  Optional integer. 1-indexed rows to keep. `NULL` keeps all.

- table:

  String. Target table, e.g. `"obs"`.

- idx_col:

  String. Index column, e.g. `"cell_idx"`.

- id_col:

  String. Identifier column, e.g. `"cell_id"`.

## Value

Invisible `NULL`. Writes the table as a side effect.
