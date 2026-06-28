# Helper to get the cache slot

Helper to get the cache slot

## Usage

``` r
.cache_slot_from_modality(modality)
```

## Arguments

- modality:

  String. One of `c("rna", "adt", "atac")`. Assumed already validated by
  the caller via `match.arg`.

## Value

The name of the cache property.
