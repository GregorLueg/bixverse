# Read how the norm of a store derives

**\[experimental\]** Reads the header only.

## Usage

``` r
rs_store_norm_kind(path)
```

## Arguments

- path:

  String. A cell- or gene-based binary (v4 or v5).

## Value

One of `"none"`, `"log_cpm"` or `"tf_idf"`. A v4 file is `"log_cpm"`.
