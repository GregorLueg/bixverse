# Check the binary format of a store before opening it

v3 stores no longer open and error with a pointer to
[`upgrade_sc_store()`](https://gregorlueg.github.io/bixverse/reference/upgrade_sc_store.md).
v4 stores still work and get a warning.

## Usage

``` r
.check_store_version(rust_con)
```

## Arguments

- rust_con:

  The `SingleCellCountData` Rust connection.

## Value

Invisible `NULL`.
