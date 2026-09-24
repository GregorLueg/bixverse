# Remove the fitted residual model

Remove the fitted residual model

## Usage

``` r
remove_residual_fit(x, ...)

## S7 method for class <bixverse::MetaCells>
remove_residual_fit(x, ...)

# S3 method for class 'ScCache'
remove_residual_fit(x, ...)

## S7 method for class <bixverse::SingleCells>
remove_residual_fit(x, ...)
```

## Arguments

- x:

  An object from which to remove the fitted model.

- ...:

  Other parameters.

## Value

The object with the fitted model dropped.

## Examples

``` r
# drop the fitted model again
sc <- demo_single_cells()
sc <- fit_residuals_sc(sc, method = "analytic_pearson", .verbose = FALSE)
sc <- remove_residual_fit(sc)
is.null(get_residual_fit(sc))
#> Warning: No fitted residual model found in the class. Run fit_residuals_sc() first. Returning NULL.
#> [1] TRUE

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
