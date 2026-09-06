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

  The class, see
  [`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md).
  Ideally, you should run
  [`preprocess_bulk_coexp()`](https://gregorlueg.github.io/bixverse/reference/preprocess_bulk_coexp.md)
  before applying this function.

- background_mat:

  Numerical matrix. The background data set.

- cor_method:

  String. Option of `c("pearson", "spearman")`.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

The class with added data to the properties for subsequent usage.

## Examples

``` r
# differential correlation of two sample groups of the same matrix
sig <- synthetic_signal_matrix()
mat <- t(sig$mat)
target <- mat[sig$group %in% c("group1", "group2"), ]
background <- mat[sig$group == "group3", ]
meta <- data.table::data.table(sample_id = rownames(target))
obj <- BulkCoExp(target, meta)
obj <- preprocess_bulk_coexp(obj, hvg = 0.3, .verbose = FALSE)
obj <- diffcor_module_processing(
  obj, background, cor_method = "pearson", .verbose = FALSE
)
obj@params$correlation_params$no_intersecting_features
#> [1] 300
```
