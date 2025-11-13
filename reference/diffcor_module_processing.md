# Prepare differential correlation-based module detection

This function will calculate the differential correlation between the
stored data set in the class and another background data set. To do so,
it uses a Fisher transformation of the correlation coefficients and
calculates a Z score based on the delta. The function will automatically
subset into shared features between the two data sets.

## Usage

``` r
diffcor_module_processing(
  object,
  background_mat,
  cor_method = c("pearson", "spearman"),
  .verbose = TRUE
)
```

## Arguments

- object:

  The class, see [`bulk_coexp()`](bulk_coexp.md). Ideally, you should
  run [`preprocess_bulk_coexp()`](preprocess_bulk_coexp.md) before
  applying this function.

- background_mat:

  Numerical matrix. The background data set.

- cor_method:

  String. Option of `c("pearson", "spearman")`.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

The class with added data to the properties for subsequent usage.
