# Update the correlation matrix to a TOM

This function will update the correlation matrix to a topological
overlap matrix. It defaults to `"v2"` and the signed version, please see
[`calculate_tom()`](https://gregorlueg.github.io/bixverse/reference/calculate_tom.md)
for details.

## Usage

``` r
cor_module_tom(object, signed = TRUE, version = c("v2", "v1"), .verbose = TRUE)
```

## Arguments

- object:

  The class, see
  [`BulkCoExp()`](https://gregorlueg.github.io/bixverse/reference/BulkCoExp.md).
  You need to have applied
  [`cor_module_processing()`](https://gregorlueg.github.io/bixverse/reference/cor_module_processing.md)
  before applying this function.

- signed:

  Boolean. Do you want to use the signed or unsigned version. Defaults
  to `TRUE`.

- version:

  String. One of `c("v2", "v1")`. Defaults to `"v2"`.

- .verbose:

  Boolean. Controls verbosity of the function.

## Value

The class with added data to the properties for subsequent usage.

## Examples

``` r
# replace the stored correlations with a signed topological overlap
mat <- t(synthetic_signal_matrix()$mat)
obj <- BulkCoExp(mat, data.table::data.table(sample_id = rownames(mat)))
obj <- preprocess_bulk_coexp(obj, hvg = 0.3, .verbose = FALSE)
obj <- cor_module_processing(obj, cor_method = "spearman", .verbose = FALSE)
obj <- cor_module_tom(obj, signed = TRUE, version = "v2", .verbose = FALSE)
obj@params$correlation_params$TOM
#> [1] TRUE
```
