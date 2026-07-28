# Get the NMF sample activity

Getter function to extract the sample activity matrix (samples x k) from
a bulk NMF fit. If no NMF fit is present, returns `NULL` with a warning.

## Usage

``` r
get_nmf_sample_activity(object)
```

## Arguments

- object:

  The class, see
  [`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md).

## Value

A samples x k numeric matrix (if found) or `NULL`.
