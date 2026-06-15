# Get the PCA loadings

Returns the PCA loadings (feature-based scores). This function is used
for the single cell-related classes and methods.

## Usage

``` r
get_pca_loadings(x)

## S7 method for class <bixverse::MetaCells>
get_pca_loadings(x)

# S3 method for class 'ScCache'
get_pca_loadings(x)

## S7 method for class <bixverse::SingleCells>
get_pca_loadings(x)
```

## Arguments

- x:

  An object to get PCA loadings from.

## Value

The PCA feature loadings from the object (if found).
