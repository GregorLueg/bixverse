# Set/add PCA factors

Set/add PCA factors

## Usage

``` r
set_pca_factors(x, pca_factor, ...)

## S7 method for class <bixverse::MetaCells>
set_pca_factors(x, pca_factor, ...)

# S3 method for class 'ScCache'
set_pca_factors(x, pca_factor, ...)

## S7 method for class <bixverse::SingleCells>
set_pca_factors(x, pca_factor, ...)

## S7 method for class <bixverse::SingleCellsSubset>
set_pca_factors(x, pca_factor, ...)
```

## Arguments

- x:

  An object to add the PCA factors for.

- pca_factor:

  Numerical matrix. The matrix with the PCA factors.

- ...:

  Other parameters.

## Examples

``` r
# an externally computed embedding pushed into the cache
sc <- demo_single_cells(prepped = FALSE)
sc <- set_pca_factors(sc, matrix(stats::rnorm(dim(sc)[1] * 2), ncol = 2))
dim(get_pca_factors(sc))
#> [1] 500   2

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
