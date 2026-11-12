# Drop a table and the named ENUM types its loads created

The named types of
[`.parquet_select_col()`](https://gregorlueg.github.io/bixverse/reference/dot-parquet_select_col.md)
are called `<table>_<column>_enum`. Only the names built from the
table's current columns are dropped, so a table whose name prefixes
another (`var`, `var_adt`) never touches the other's types. A type
cannot be dropped while a table uses it, hence the table goes first.

## Usage

``` r
.drop_table_and_enums(con, table)
```

## Arguments

- con:

  DuckDB connection.

- table:

  String. The table, e.g. `"obs"`.

## Value

Invisible `NULL`.
