# Write the cache for a given modality back onto the object

Write the cache for a given modality back onto the object

## Usage

``` r
.sc_cache_put(x, cache, modality = "rna")
```

## Arguments

- x:

  The object owning the cache.

- cache:

  The cache list to write.

- modality:

  String. One of `c("rna", "adt", "atac", "wnn")`.

## Value

The object with the cache written back.
