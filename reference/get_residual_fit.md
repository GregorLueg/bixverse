# Get the fitted residual model

Returns the `ScResidualFit` written by
[`fit_residuals_sc()`](https://gregorlueg.github.io/bixverse/reference/fit_residuals_sc.md).
This function is used for the single cell-related classes and methods.

Warns and returns `NULL` when nothing was fitted, so it is safe to use
as a presence probe. The functions that compute with the fit assert
instead.

## Usage

``` r
get_residual_fit(x, ...)

## S7 method for class <bixverse::MetaCells>
get_residual_fit(x, ...)

# S3 method for class 'ScCache'
get_residual_fit(x, ...)

## S7 method for class <bixverse::SingleCells>
get_residual_fit(x, ...)
```

## Arguments

- x:

  An object to get the fitted model from.

- ...:

  Other parameters.

## Value

The `ScResidualFit` object, or `NULL` when nothing was fitted.

## Examples

``` r
# the fitted model and the genes it covers
sc <- demo_single_cells()
sc <- fit_residuals_sc(sc, method = "analytic_pearson", .verbose = FALSE)
length(get_residual_fit(sc)$genes)
#> [1] 50

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
