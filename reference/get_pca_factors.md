# Get the PCA factors

Returns the PCA factors (sample-based scores). This function is used for
the single cell-related classes and methods.

## Usage

``` r
get_pca_factors(x, ...)

## S7 method for class <bixverse::MetaCells>
get_pca_factors(x, ...)

# S3 method for class 'ScCache'
get_pca_factors(x, ...)

## S7 method for class <bixverse::SingleCells>
get_pca_factors(x, ...)
```

## Arguments

- x:

  An object to get PCA factors from.

- ...:

  Other parameters.

## Value

The PCA factors from the object (if found).
