# Largest categorical cast to an inline ENUM literal

Up to this many levels, a categorical is cast with `ENUM('a', 'b', ...)`
written into the SQL, which is what
[`DBI::dbWriteTable()`](https://dbi.r-dbi.org/reference/dbWriteTable.html)
produces for a factor. Beyond it, parsing the literal dominates the load
(a per-cell categorical can carry millions of levels), so the levels go
in through a registered data.frame and a named type instead. Not tuned.

## Usage

``` r
SC_ENUM_LITERAL_MAX
```
