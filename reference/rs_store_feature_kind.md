# Read what the features of a store are

**\[experimental\]** Reads the header only.

## Usage

``` r
rs_store_feature_kind(path)
```

## Arguments

- path:

  String. A cell- or gene-based binary (v4 or v5).

## Value

One of `"gene"`, `"peak"`, `"tile"`, `"snp"`, `"protein"` or
`"gene_activity"`. A v4 file is `"gene"`.
