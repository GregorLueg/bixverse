# Set/add the fitted residual model

Set/add the fitted residual model

## Usage

``` r
set_residual_fit(x, residual_fit, ...)

## S7 method for class <bixverse::MetaCells>
set_residual_fit(x, residual_fit, ...)

# S3 method for class 'ScCache'
set_residual_fit(x, residual_fit, ...)

## S7 method for class <bixverse::SingleCells>
set_residual_fit(x, residual_fit, ...)
```

## Arguments

- x:

  An object to add the fitted model to.

- residual_fit:

  `ScResidualFit` class, as returned by
  [`fit_residuals_sc()`](https://gregorlueg.github.io/bixverse/reference/fit_residuals_sc.md).

- ...:

  Other parameters.

## Value

The object with the fitted model attached.

## Examples

``` r
# the fitted model taken out and put back
sc <- demo_single_cells()
sc <- fit_residuals_sc(sc, method = "analytic_pearson", .verbose = FALSE)
fit <- get_residual_fit(sc)
sc <- set_residual_fit(remove_residual_fit(sc), fit)
get_residual_fit(sc)$method
#> [1] "analytic_pearson"

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
