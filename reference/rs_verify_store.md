# Check every chunk of a store file against its checksums

**\[experimental\]**

## Usage

``` r
rs_verify_store(path)
```

## Arguments

- path:

  String. A cell- or gene-based binary.

## Value

`TRUE` when every chunk matched, `FALSE` for a file without checksums
(v4). Errors on a corrupt chunk.
