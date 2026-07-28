# Get the NMF module membership data.table

Getter function to extract the gene-to-module data.table from a bulk NMF
fit. Each row is one gene assigned to its top-loading module.

## Usage

``` r
get_nmf_modules(object)
```

## Arguments

- object:

  The class, see
  [`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md).

## Value

A data.table with `gene`, `module_id`, `loading`, `sign` columns (if
found) or `NULL`.
