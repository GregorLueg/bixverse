# The cell-ID hash in a binary's header

**\[experimental\]**

## Usage

``` r
rs_cell_id_hash(path)
```

## Arguments

- path:

  String. A cell- or gene-based binary.

## Value

The 16-byte hash as hex, `NULL` when the store records no cell IDs.
