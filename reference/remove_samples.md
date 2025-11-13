# Remove samples from object

This function allows to remove certain samples from the object

## Usage

``` r
remove_samples(object, samples_to_remove, ...)
```

## Arguments

- object:

  The underlying object, either [`bixverse::bulk_coexp`](bulk_coexp.md)
  or [`bixverse::bulk_dge`](bulk_dge.md).

- samples_to_remove:

  Character vector. The sample identifiers to remove.

- ...:

  Additional arguments to parse to the functions.

## Value

Returns the object with the samples removed. This will regenerated the
object from the start and remove any data in it.
