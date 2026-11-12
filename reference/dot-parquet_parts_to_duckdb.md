# Load several Parquet metadata tables into one DuckDB table

The multi-file counterpart of
[`.parquet_to_duckdb()`](https://gregorlueg.github.io/bixverse/reference/dot-parquet_to_duckdb.md).
Stacks the parts in order with `UNION ALL BY NAME` (columns missing from
a part become `NULL`), adds `exp_id`, prefixes the identifier with it
and renumbers the index. Column order is the index, the identifier,
`exp_id`, then the remaining columns in order of first appearance, as
`data.table::rbindlist(fill = TRUE)` produced.

## Usage

``` r
.parquet_parts_to_duckdb(
  con,
  table,
  parts,
  idx_col,
  id_col,
  enum_literal_max = SC_ENUM_LITERAL_MAX
)
```

## Arguments

- con:

  DuckDB connection.

- table:

  String. Target table, e.g. `"obs"`.

- parts:

  List. One element per file, each a list with `path` (the Parquet
  file), `info` (its writer output) and `exp_id`.

- idx_col:

  String. The index column, e.g. `"cell_idx"`.

- id_col:

  String. The identifier column, e.g. `"cell_id"`.

- enum_literal_max:

  Integer. Most levels cast with an inline ENUM, see
  [SC_ENUM_LITERAL_MAX](https://gregorlueg.github.io/bixverse/reference/SC_ENUM_LITERAL_MAX.md).

## Value

Invisible `NULL`. Writes the table as a side effect.
