# Remove samples from object

This function allows to remove certain samples from the object

## Usage

``` r
remove_samples(object, samples_to_remove, ...)
```

## Arguments

- object:

  The underlying object, either
  [`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md)
  or
  [`BulkDge()`](https://gregorlueg.github.io/bixverse/reference/BulkDge.md).

- samples_to_remove:

  Character vector. The sample identifiers to remove.

- ...:

  Additional arguments to parse to the functions.

## Value

Returns the object with the samples removed. This will regenerated the
object from the start and remove any data in it.

## Examples

``` r
# drop the first two samples and rebuild the object
syn <- synthetic_bulk_cor_matrix()
meta <- data.table::data.table(sample_id = colnames(syn$counts))
object <- BulkDge(raw_counts = syn$counts, meta_data = meta)
object <- remove_samples(object, c("sample_1", "sample_2"))
dim(S7::prop(object, "raw_counts"))
#> [1] 1000   98
```
