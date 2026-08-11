# Read the cache for a given modality

Read the cache for a given modality

## Usage

``` r
.sc_cache_get(x, modality = "rna")
```

## Arguments

- x:

  The object owning the cache.

- modality:

  String. One of `c("rna", "adt", "atac", "wnn")`.

## Value

The cache list, or `NULL` when the object has no such cache.
