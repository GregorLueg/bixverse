# SQL select expression for one Parquet column

A categorical with at least one level is cast to an ENUM with its levels
in order, which is what
[`DBI::dbWriteTable()`](https://dbi.r-dbi.org/reference/dbWriteTable.html)
makes of a factor. `ENUM()` with no levels is invalid, so an empty
categorical stays VARCHAR. Up to `enum_literal_max` levels the ENUM is
written inline; above it a named type `<table>_<column>_enum` is created
from a registered data.frame, which is far cheaper for very many levels.
Either way the column comes back into R as a factor with the same
levels.

## Usage

``` r
.parquet_select_col(
  con,
  table,
  name,
  levels,
  enum_literal_max = SC_ENUM_LITERAL_MAX
)
```

## Arguments

- con:

  DuckDB connection.

- table:

  String. The table being built, used to name the type.

- name:

  String. Column name.

- levels:

  Character vector or `NULL`. The ENUM levels.

- enum_literal_max:

  Integer. Most levels written inline.

## Value

String. The select expression. Creates the named type as a side effect
where one is needed.
