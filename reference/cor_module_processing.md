# Prepare correlation-based module detection

This function will calculate the correlation coefficients between the
genes, using the highly variable genes (if available, otherwise the
function will use the raw data). The data will be stored in a
memory-efficient format in the properties of the class.

## Usage

``` r
cor_module_processing(
  object,
  cor_method = c("pearson", "spearman"),
  .verbose = TRUE
)
```

## Arguments

- object:

  The class, see [`bulk_coexp()`](bulk_coexp.md). Ideally, you should
  run [`preprocess_bulk_coexp()`](preprocess_bulk_coexp.md) before
  applying this function.

- cor_method:

  String. Option of `c("pearson", "spearman")`.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

The class with added data to the properties for subsequent usage.
