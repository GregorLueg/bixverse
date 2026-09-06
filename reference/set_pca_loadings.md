# Set/add PCA loadings

Set/add PCA loadings

## Usage

``` r
set_pca_loadings(x, pca_loading, ...)

## S7 method for class <bixverse::MetaCells>
set_pca_loadings(x, pca_loading, ...)

# S3 method for class 'ScCache'
set_pca_loadings(x, pca_loading, ...)

## S7 method for class <bixverse::SingleCells>
set_pca_loadings(x, pca_loading, ...)

## S7 method for class <bixverse::SingleCellsSubset>
set_pca_loadings(x, pca_loading, ...)
```

## Arguments

- x:

  An object to add the PCA loadings for.

- pca_loading:

  Numerical matrix. The Matrix with the PCA loadings.

- ...:

  Other parameters.

## Examples

``` r
# only the first two loading vectors kept
sc <- demo_single_cells()
sc <- set_pca_loadings(sc, get_pca_loadings(sc)[, 1:2])
dim(get_pca_loadings(sc))
#> [1] 30  2

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
