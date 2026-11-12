# Load a Parquet metadata table into DuckDB

Runs `CREATE OR REPLACE TABLE ... AS SELECT ... FROM read_parquet(...)`,
casting the categorical columns to ENUMs.

## Usage

``` r
.parquet_to_duckdb(
  con,
  table,
  parquet_path,
  info,
  enum_literal_max = SC_ENUM_LITERAL_MAX
)
```

## Arguments

- con:

  DuckDB connection.

- table:

  String. Target table, e.g. `"obs"`.

- parquet_path:

  String. The Parquet file.

- info:

  List. Output of
  [`rs_h5ad_frame_to_parquet()`](https://gregorlueg.github.io/bixverse/reference/rs_h5ad_frame_to_parquet.md)
  or
  [`rs_tenx_frame_to_parquet()`](https://gregorlueg.github.io/bixverse/reference/rs_tenx_frame_to_parquet.md)
  for that file.

- enum_literal_max:

  Integer. Most levels cast with an inline ENUM, see
  [SC_ENUM_LITERAL_MAX](https://gregorlueg.github.io/bixverse/reference/SC_ENUM_LITERAL_MAX.md).

## Value

Invisible `NULL`. Writes the table as a side effect.
