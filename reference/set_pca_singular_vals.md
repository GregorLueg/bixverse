# Set/add PCA singular values

Set/add PCA singular values

## Usage

``` r
set_pca_singular_vals(x, singular_vals, ...)

## S7 method for class <bixverse::MetaCells>
set_pca_singular_vals(x, singular_vals, ...)

# S3 method for class 'ScCache'
set_pca_singular_vals(x, singular_vals, ...)

## S7 method for class <bixverse::SingleCells>
set_pca_singular_vals(x, singular_vals, ...)

## S7 method for class <bixverse::SingleCellsSubset>
set_pca_singular_vals(x, singular_vals, ...)
```

## Arguments

- x:

  An object to add the singular values for.

- singular_vals:

  Numerical vector. The singular values.

- ...:

  Other parameters.

## Examples

``` r
# singular values from a decomposition done elsewhere
sc <- demo_single_cells(prepped = FALSE)
sc <- set_pca_singular_vals(sc, c(4.1, 2.3))
get_pca_singular_val(sc)
#> [1] 4.1 2.3

unlink(sc@dir_data, recursive = TRUE, force = TRUE)
```
