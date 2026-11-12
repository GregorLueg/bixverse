# ENUM levels of a column across several Parquet parts

Reproduces what
[`data.table::rbindlist()`](https://rdrr.io/pkg/data.table/man/rbindlist.html)
makes of the column when at least one part holds it as a factor: per
part in order, a categorical adds its levels in level order and a string
column adds its values in row order, each only when not yet seen.

## Usage

``` r
.parquet_parts_levels(con, name, parts, idx_col)
```

## Arguments

- con:

  DuckDB connection.

- name:

  String. Column name.

- parts:

  List of parts, see
  [`.parquet_parts_to_duckdb()`](https://gregorlueg.github.io/bixverse/reference/dot-parquet_parts_to_duckdb.md).

- idx_col:

  String. The per-part row index column.

## Value

Character vector of levels, or `NULL` if the column is not categorical
in any part, or is held as something other than a string in some part.
